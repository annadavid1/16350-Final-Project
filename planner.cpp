#include <math.h>
#include <iostream>
#include <SFML/Graphics.hpp>
#include <thread>
#include <mutex>
#include <cmath>
#include <random>
#include <queue>
#include <chrono>
#include <vector>

using namespace std;

#define G 9.8
#define PI M_PI
#define DOUBLEMAX std::numeric_limits<double>::max()

#define WINDOWX 800
#define WINDOWY 600
#define WINDOWSCALE 100

#define VISX(x) WINDOWSCALE*x + WINDOWX/2
#define VISY(y) WINDOWY - WINDOWSCALE*y - 50
#define BALLRADIUS 0.2

#define BASEX 0.0
#define BASEY 0.0
#define LINK1 1.5
#define LINK2 1.5

#define EFFECTORRADIUS 0.1
#define CONTACTRADIUS EFFECTORRADIUS + BALLRADIUS

#define TARGETX 0.0
#define TARGETY BASEY + LINK1
#define TARGETRADIUS (LINK2 * 3.0 / 4.0)
#define MINHEIGHT BASEY + LINK1 + 2.0

#define MAXOMEGA1 2*PI
#define MAXOMEGA2 2*PI
#define MAXACC1 4*PI
#define MAXACC2 4*PI

// epsilon for extend in RRT
#define EPSILON PI/8.0

#define SAMPLETIME 2.25

#define NUMSAMPLES 32

// extend outputs
#define TRAPPED 0
#define ADVANCED 1
#define REACHED 2


// initialize random distribution
random_device rd;
default_random_engine gen(rd());
uniform_real_distribution<> rand_omega1(-MAXOMEGA1, MAXOMEGA1);
uniform_real_distribution<> rand_omega2(-MAXOMEGA2, MAXOMEGA2);
bernoulli_distribution sampleBiased(0.7);


struct BallState
{
    double x;
    double y;
    double vx;
    double vy;
    double t;
};

struct CatchTarget
{
    double x;
    double y;
    double t;
};

struct OmegaState
{
    double omega1;
    double omega2;
    double t;
};

struct AngleState
{
    double theta1;
    double theta2;
    double omega1;
    double omega2;
    double t;
};

struct EffectorState
{
    double xmid;
    double ymid;
    double x;
    double y;
    double vx;
    double vy;
    double t;
};

struct WorldState
{
    AngleState angleStart;
    AngleState angleEnd;
    EffectorState effector;
    BallState ballStart;
    BallState ballNext;
    float time;
};

struct Node
{
    AngleState angles;
    Node *bp;

    Node(AngleState angles) : angles(angles), bp(nullptr) {}
    Node(AngleState angles, Node *bp) : angles(angles), bp(bp) {}
};

// tree storing root and list of nodes
struct Tree {
	Node *root; // start pos
	vector<Node*> V; // list of nodes in tree for easy iterating

	// tree initialized with just root
	Tree(Node *root) : root(root), V({root}) {};
};

WorldState world;
mutex worldLock;

queue<pair<Node*, BallState>> planQueue;
mutex queueLock;

double startTime = std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::system_clock::now().time_since_epoch()).count() / 1000.0;


// get position of ball at time t that started at start
static BallState ball(BallState start, double t) {
    double dt = t - start.t;
    double x = start.x + start.vx * dt;
    double y = start.y + start.vy * dt - G * dt * dt / 2;
    double vx = start.vx;
    double vy = start.vy - G * dt;
    return {x, y, vx, vy, t};
}

// get position of arm at time t moving from start to end
static AngleState arm(AngleState start, AngleState end, double t) {
    double tf = end.t - start.t;
    double dt = t - start.t;
    if (dt < 0 || tf <= 0) {
        return start;
    }
    double a1 = (end.omega1 - start.omega1)/tf;
    double omega1 = start.omega1 + a1*dt;
    double theta1 = start.theta1 + start.omega1*dt + a1*dt*dt/2;
    double a2 = (end.omega2 - start.omega2)/tf;
    double omega2 = start.omega2 + a2*dt;
    double theta2 = start.theta2 + start.omega2*dt + a2*dt*dt/2;
    return {theta1, theta2, omega1, omega2, t};
}

// get position of arm at time t moving from start to end
static AngleState arm(AngleState start, OmegaState end, double t) {
    double tf = end.t - start.t;
    double dt = t - start.t;
    if (dt < 0 || tf <= 0) {
        return start;
    }
    double a1 = (end.omega1 - start.omega1)/tf;
    double omega1 = start.omega1 + a1*dt;
    double theta1 = start.theta1 + start.omega1*dt + a1*dt*dt/2;
    double a2 = (end.omega2 - start.omega2)/tf;
    double omega2 = start.omega2 + a2*dt;
    double theta2 = start.theta2 + start.omega2*dt + a2*dt*dt/2;
    return {theta1, theta2, omega1, omega2, t};
}

// get end effector state from angle state
static EffectorState anglesToEffector(AngleState angles) {
    double xmid = LINK1 * cos(angles.theta1);
    double x = xmid + LINK2 * cos(angles.theta1+angles.theta2);
    double ymid = LINK1 * sin(angles.theta1);
    double y = ymid + LINK2 * sin(angles.theta1+angles.theta2);
    double vx = -LINK1 * sin(angles.theta1) * angles.omega1 
        - LINK2 * sin(angles.theta1 + angles.theta2) * (angles.omega1 + angles.omega2);
    double vy = LINK1 * cos(angles.theta1) * angles.omega1
        + LINK2 * cos(angles.theta1 + angles.theta2) * (angles.omega1 + angles.omega2);
    return {xmid, ymid, x, y, vx, vy, angles.t};
}

static pair<bool, BallState> isGoalConfig(EffectorState effector, BallState ballStart) {
    if (effector.vy <= 0) {
        return {false, {0,0,0,0,0}};
    }
    BallState currBall = ball(ballStart, effector.t);
    int i;
    double t, t0, xdist, ydist, dist, x0, y0, vx, vy, a, b, c, d, e, A, B, C, D,
        p, q, u, v, r, phi, m, minT, maxDist;
    double ts[3];
    xdist = currBall.x - effector.x;
    ydist = currBall.y - effector.y;
    dist = xdist * xdist + ydist * ydist;
    maxDist = EFFECTORRADIUS + BALLRADIUS;
    if (dist > maxDist * maxDist) {
        // cout << "not goal 1\n";
        // cout << ydist << "," << xdist << "," << currBall.y << "," << effector.y << "," << effector.t << endl;
        return {false, {0,0,0,0,0}};
    }
    x0 = currBall.x - TARGETX;
    y0 = currBall.y - TARGETY;
    vx = effector.vx;
    vy = effector.vy;
    t0 = effector.t;

    if (y0 + vy*vy/(2*G) < MINHEIGHT - TARGETY) {
        // cout << "not goal 2\n";
        return {false, {0,0,0,0,0}};
    }
    minT = vy / G;
    
    // want at^4 + bt^3 + ct^2 + dt + e <= 0
    // solve by finding roots of derivative: 4at^3 + 3bt^2 + 2ct + d = 0
    a = G*G/4.0;
    b = -vy*G;
    c = vx*vx + vy*vy - y0*G;
    d = 2.0*x0*vx + 2.0*y0*vy;
    e = x0*x0 + y0*y0 - TARGETRADIUS*TARGETRADIUS;

    // get into t^3 + At^2 + Bt + C = 0
    A = 3.0*b / (4.0*a);
    B = c / (2.0*a);
    C = d / (4.0*a);

    // use t = s - A/3 to get into s^3 + ps + q = 0
    p = B - A*A/3.0;
    q = 2*A*A*A/27.0 - A*B/3.0 + C;

    // discriminant
    D = q*q/4.0 + p*p*p/27.0;

    if (D >= 0) {
        // one real root
        u = cbrt(-q/2.0 + sqrt(D));
        v = cbrt(-q/2.0 - sqrt(D));
        t = u + v - A/3.0;
        if (t < minT) {
            // cout << "not goal 3\n";
            return {false, {0,0,0,0,0}};
        }
        xdist = x0+vx*t;
        ydist = y0+vy*t-G*t*t/2;
        if (xdist*xdist + ydist*ydist <= TARGETRADIUS*TARGETRADIUS) {
            // cout << "definitely goal 4\n";
            return {true, {x0 + TARGETX, y0 + TARGETY, vx, vy, t0}};
        }
        // cout << "not goal 5\n";
        return {false, {0,0,0,0,0}};
    }

    r = sqrt(-p*p*p/27.0);
    phi = acos(-q/(2.0*r));
    m = 2 * sqrt(-p/3.0);

    ts[0] = m * cos(phi/3.0) - A/3.0;
    ts[1] = m * cos((phi + 2.0*PI)/3.0) - A/3.0;
    ts[2] = m * cos((phi + 4.0*PI)/3.0) - A/3.0;
    for (i = 0; i < 3; i++) {
        t = ts[i];
        if (t < minT) {
            continue;
        }
        double res = a*t*t*t*t + b*t*t*t + c*t*t + d*t + e;
        if (res <= 0) {
            // cout << "definitely goal 6\n";
            return {true, {x0 + TARGETX, y0 + TARGETY, vx, vy, t0}};
        }
    }
    // cout << "not goal 7\n";
    return {false, {0,0,0,0,0}};
}

// check whether a given effector state can successfully catch/throw ball
static pair<bool, BallState> isGoalConfig(AngleState angles, BallState ballStart) {
    return isGoalConfig(anglesToEffector(angles), ballStart);
}

static pair<bool, BallState> isGoalConfig2(AngleState angles, BallState ballStart, AngleState& anglesp) {
    double T_MAX = 1, DT = 0.1, DA = MAXACC1 / 16, theta1p, theta2p, omega1p, omega2p, dx, dy;
    EffectorState effectorp;
    BallState ballt;
    pair<bool, BallState> isGoal;
    for (double t = 0; t <= T_MAX; t += DT) {
        ballt = ball(ballStart, t+angles.t);
        for (double a1 = -MAXACC1; a1 <= MAXACC1; a1 += DA) {
            for (double a2 = -MAXACC2; a2 <= MAXACC2; a2 += DA) {
                // compute state
                theta1p = angles.theta1 + angles.omega1*t + 0.5*a1*t*t;
                theta2p = angles.theta2 + angles.omega2*t + 0.5*a2*t*t;
                omega1p = angles.omega1 + a1*t;
                omega2p = angles.omega2 + a2*t;
                anglesp = {theta1p, theta2p, omega1p, omega2p, t + angles.t};

                isGoal = isGoalConfig(anglesp, ballStart);
                if (isGoal.first) {
                    return isGoal;
                }

                // effectorp = anglesToEffector(anglesp);

                // // constraint 4
                // dx = effectorp.x - ballt.x;
                // dy = effectorp.y - ballt.y;
                // if (dx*dx + dy*dy > CONTACTRADIUS*CONTACTRADIUS) {continue;}

                // // constraint 6
                // if (effectorp.y + effectorp.vy*effectorp.vy/(2*G) < MINHEIGHT) {continue;}

                // // constraint 5
                // if (hitsTarget(...)) {
                //     return true;
                // }
            }
        }
    }
    return {false, {0,0,0,0,0}};
}

// Euclidean distance between 2 angle configs
static double dist(AngleState start, AngleState end) {
	double omega1 = start.omega1 - end.omega1;
    double omega2 = start.omega2 - end.omega2;
    double t = (start.t - end.t)/4; // scale time down so it doesn't overpower
	return sqrt(omega1*omega1 + omega2*omega2 + t*t);
}

// Euclidean distance between angle and omega configs
static double dist(AngleState start, OmegaState end) {
	double omega1 = start.omega1 - end.omega1;
    double omega2 = start.omega2 - end.omega2;
    double t = (start.t - end.t)/4; // scale time down so it doesn't overpower
	return sqrt(omega1*omega1 + omega2*omega2 + t*t);
}

double distToBall(AngleState q, BallState ballStart) {
    BallState b = ball(ballStart, q.t);
    EffectorState e = anglesToEffector(q);
    double dx = e.x - b.x;
    double dy = e.y - b.y;
    return sqrt(dx*dx + dy*dy);
}

double distToBall(AngleState start, OmegaState q, BallState ballStart) {
    return distToBall(arm(start, q, q.t), ballStart);
}

// check end after start and acceleration less than max
static bool canConnect(AngleState start, AngleState end) {
    double t = end.t - start.t;
    if (t <= 0) {
        return false;
    }
    double a1 = (end.omega1 - start.omega1)/t;
    double a2 = (end.omega2 - start.omega2)/t;
    return -MAXACC1 <= a1 && a1 <= MAXACC1 && -MAXACC2 <= a2 && a2 <= MAXACC2;
}

// check end after start and acceleration less than max
static bool canConnect(AngleState start, OmegaState end) {
    double t = end.t - start.t;
    if (t <= 0) {
        return false;
    }
    double a1 = (end.omega1 - start.omega1)/t;
    double a2 = (end.omega2 - start.omega2)/t;
    return -MAXACC1 <= a1 && a1 <= MAXACC1 && -MAXACC2 <= a2 && a2 <= MAXACC2;
}

// get random angular velocities and time after startT
OmegaState sample(double startT, double endT) {
    double omega1 = rand_omega1(gen);
    double omega2 = rand_omega2(gen);
    uniform_real_distribution<> rand_time(startT, endT);
    double t = rand_time(gen);
    return {omega1, omega2, t};
}

CatchTarget sampleCatchTarget(BallState ballStart) {
    uniform_real_distribution<> rand_time(
        ballStart.t,
        ballStart.t + ballStart.vy / G * SAMPLETIME
    );
    double t = rand_time(gen);
    BallState b = ball(ballStart, t);
    return {b.x, b.y, t};
}

double distToTarget(AngleState q, CatchTarget target) {
    EffectorState e = anglesToEffector(q);
    double dx = e.x - target.x;
    double dy = e.y - target.y;
    double dt = (q.t - target.t) / 4.0;
    return sqrt(dx * dx + dy * dy + dt * dt);
}

bool solveIK(double x, double y, double& theta1, double& theta2, double& theta12, double& theta22) {
    double r2 = x * x + y * y;
    double c2 = (r2 - LINK1 * LINK1 - LINK2 * LINK2) / (2 * LINK1 * LINK2);
    if (c2 < -1.0 || c2 > 1.0) return false;

    double s2 = sqrt(1.0 - c2 * c2); // elbow-down branch
    theta2 = atan2(s2, c2);

    double k1 = LINK1 + LINK2 * c2;
    double k2 = LINK2 * s2;
    theta1 = atan2(y, x) - atan2(k2, k1);

    theta22 = atan2(-s2, c2);

    k1 = LINK1 + LINK2 * c2;
    k2 = -LINK2 * s2;
    theta12 = atan2(y, x) - atan2(k2, k1);
    return true;
}

bool makeGoalState(
    AngleState start,
    CatchTarget target,
    AngleState& goal
) {
    double theta1, theta2, theta12, theta22;
    AngleState goal2;
    if (!solveIK(target.x, target.y, theta1, theta2, theta12, theta22)) return false;

    double dt = target.t - start.t;
    if (dt <= 0) return false;

    double omega1 = 2 * (theta1 - start.theta1) / dt - start.omega1;
    double omega2 = 2 * (theta2 - start.theta2) / dt - start.omega2;
    
    goal = {theta1, theta2, omega1, omega2, target.t};

    double omega12 = 2 * (theta12 - start.theta1) / dt - start.omega1;
    double omega22 = 2 * (theta22 - start.theta2) / dt - start.omega2;
    
    goal2 = {theta12, theta22, omega12, omega22, target.t};

    if (anglesToEffector(goal2).vy > anglesToEffector(goal).vy) {
        AngleState temp = goal;
        goal = goal2;
        goal2 = temp;
        double temp2 = omega1;
        omega1 = omega12;
        omega12 = temp2;
        temp2 = omega2;
        omega2 = omega22;
        omega22 = temp2;
    }

    if (omega1 > -MAXOMEGA1 && MAXOMEGA1 > omega1 && omega2 > -MAXOMEGA2 
        && MAXOMEGA2 > omega2 && canConnect(start, goal)) {
        return true;
    }
    if (omega12 > -MAXOMEGA1 && MAXOMEGA1 > omega12 && omega22 > -MAXOMEGA2 
        && MAXOMEGA2 > omega22 && canConnect(start, goal2)) {
        AngleState temp = goal;
        goal = goal2;
        goal2 = temp;
        double temp2 = omega1;
        omega1 = omega12;
        omega12 = temp2;
        temp2 = omega2;
        omega2 = omega22;
        omega22 = temp2;
        return true;
    }

    return false;
}

Node* nearestToTarget(Tree* t, CatchTarget target) {
    Node* qmin = nullptr;
    double minDist = DOUBLEMAX;
    AngleState goal;
    for (Node* q : t->V) {
        double d = distToTarget(q->angles, target);
        if (d < minDist && makeGoalState(q->angles, target, goal)) {
            qmin = q;
            minDist = d;
        }
    }
    return qmin;
}

static pair<int, Node*> extend(Tree* t, BallState ballStart) {
    // if (sampleBiased(gen)) {
    //     CatchTarget target = sampleCatchTarget(ballStart);
    //     Node* qmin = nearestToTarget(t, target);
    //     if (qmin == nullptr) return {TRAPPED, nullptr};

    //     AngleState goal;
    //     if (!makeGoalState(qmin->angles, target, goal)) return {TRAPPED, nullptr};
    //     if (!canConnect(qmin->angles, goal)) return {TRAPPED, nullptr};

    //     double d = dist(qmin->angles, goal);
    //     double scale = EPSILON / d;

    //     AngleState qextAngle;
    //     if (scale >= 1.0) {
    //         qextAngle = goal;
    //     } else {
    //         OmegaState qext = {
    //             qmin->angles.omega1 + (goal.omega1 - qmin->angles.omega1) * scale,
    //             qmin->angles.omega2 + (goal.omega2 - qmin->angles.omega2) * scale,
    //             qmin->angles.t + (goal.t - qmin->angles.t) * scale
    //         };
    //         qextAngle = arm(qmin->angles, qext, qext.t);
    //     }

    //     Node* qnew = new Node(qextAngle, qmin);
    //     AngleState qangles, qcheck;
    //     qangles = qnew->angles;
    //     for (int i = 1; i < NUMSAMPLES; i++) {
    //         qcheck = arm(qmin->angles, qangles, qmin->angles.t 
    //             + (qangles.t - qmin->angles.t) * i / NUMSAMPLES);
    //         if (isGoalConfig(qcheck, ballStart).first) {
    //             qnew = new Node(qcheck, qmin);
    //             t->V.push_back(qnew);
    //             return {ADVANCED, qnew};
    //         }
    //     }
    //     t->V.push_back(qnew);
    //     return {REACHED, qnew};
    //     t->V.push_back(qnew);
    //     return {scale >= 1.0 ? REACHED : ADVANCED, qnew};
    // } else {
    OmegaState qrand = sample(ballStart.t, 
        ballStart.t + ballStart.vy / G * SAMPLETIME);
    Node *q, *qnew, *qmin = nullptr;
    double d;
    double minDist = DOUBLEMAX;
    for (int i = 0; i < t->V.size(); i++) {
        q = t->V[i];
        d = dist(q->angles, qrand);
        if (d < minDist && canConnect(q->angles, qrand)) {
            qmin = q;
            minDist = d;
        }
    }
    if (qmin == nullptr) {
        return {TRAPPED, nullptr};
    }
    double scale = EPSILON / minDist;
    AngleState qcheck, qangles;
    if (scale >= 1) {
        qnew = new Node(arm(qmin->angles, qrand, qrand.t), qmin);
        // qangles = qnew->angles;
        // for (int i = 1; i < NUMSAMPLES; i++) {
        //     qcheck = arm(qmin->angles, qangles, qmin->angles.t 
        //         + (qangles.t - qmin->angles.t) * i / NUMSAMPLES);
        //     if (isGoalConfig(qcheck, ballStart).first) {
        //         qnew = new Node(qcheck, qmin);
        //         t->V.push_back(qnew);
        //         return {ADVANCED, qnew};
        //     }
        // }
        t->V.push_back(qnew);
        return {REACHED, qnew};
    }

    OmegaState qext = {
        qmin->angles.omega1 + (qrand.omega1 - qmin->angles.omega1) * scale,
        qmin->angles.omega2 + (qrand.omega2 - qmin->angles.omega2) * scale,
        qmin->angles.t + (qrand.t - qmin->angles.t) * scale
    };
    qnew = new Node(arm(qmin->angles, qext, qext.t), qmin);
    qangles = qnew->angles;
    // for (int i = 1; i < NUMSAMPLES; i++) {
    //     qcheck = arm(qmin->angles, qangles, qmin->angles.t 
    //         + (qangles.t - qmin->angles.t) * i / NUMSAMPLES);
    //     if (isGoalConfig(qcheck, ballStart).first) {
    //         qnew = new Node(qcheck, qmin);
    //         t->V.push_back(qnew);
    //         return {ADVANCED, qnew};
    //     }
    // }
    t->V.push_back(qnew);
    return {ADVANCED, qnew};
    // }
}

// static pair<int, Node*> extend(Tree *t, OmegaState qrand, BallState ballStart) {
//     int n = t->V.size();
//     OmegaState qext;
//     AngleState qextAngle;
//     Node *q, *qnew, *qmin = nullptr;
//     double d, scale, minDist = DOUBLEMAX;
//     for (int i = 0; i < n; i++) {
//         q = t->V[i];
//         d = distToBall(q->angles, qrand, ballStart);
//         if (d < minDist && canConnect(q->angles, qrand)) {
//             qmin = q;
//             minDist = d;
//         }
//     }
//     if (qmin == nullptr) {
//         return {TRAPPED, nullptr};
//     }
//     scale = EPSILON / dist(qmin->angles, qrand);
//     if (scale >= 1) {
//         qextAngle = arm(qmin->angles, qrand, qrand.t);
//         qnew = new Node(qextAngle, qmin);
//         t->V.push_back(qnew);
//         return {2, qnew};
//     }
//     qext = {qmin->angles.omega1 + (qrand.omega1 - qmin->angles.omega1)*scale,
//             qmin->angles.omega2 + (qrand.omega2 - qmin->angles.omega2)*scale,
//             qmin->angles.t + (qrand.t - qmin->angles.t)*scale};

//     qextAngle = arm(qmin->angles, qext, qext.t);
//     qnew = new Node(qextAngle, qmin);
//     t->V.push_back(qnew);
//     return {1, qnew};
// }

static pair<Node*, BallState> runRRT(
    AngleState startAngle,
    BallState ballStart
) {
	int s, pts, i, j;
	double d, diff, cost = 0.0;
	OmegaState qrand, pos;
	Node *start, *qnew, *q1, *q2;
	Tree *t;
	vector<double*> path = {};
	pair<int, Node*> ext;
    pair<bool, BallState> isGoal;
    AngleState anglesp;

	// initialize start and goal nodes and trees
	start = new Node(startAngle);
	t = new Tree(start);

	// main loop, goes until any path found
	while (1) {
		// get random config and extend t towards it
		// qrand = sample(ballStart.t, ballStart.t + ballStart.vy/G*SAMPLETIME);
		ext = extend(t, ballStart);
		s = ext.first; qnew = ext.second;
		if (s != TRAPPED) {
            isGoal = isGoalConfig2(qnew->angles, ballStart, anglesp);
            if (isGoal.first) {
                // reached goal
                qnew = new Node(anglesp, qnew);
                return {qnew, isGoal.second};
            }
		}
	}
    // should never get here
}

/*
* planner thread creates plan from current (slowed) state to the first goal 
* state found, then updates the current plan, then starts planning for the 
* next goal during execution
*/
void plannerThread() {
    pair<Node*, BallState> plan;
    WorldState snapshot;
    {
        lock_guard<mutex> lock(worldLock);
        snapshot = world;
    }
    AngleState currAngles = snapshot.angleEnd;
    BallState ballStart = snapshot.ballStart;
    while (1)
    {
        cout << "Starting RRT...\n";
        plan = runRRT(currAngles, ballStart);
        {
            lock_guard<mutex> lock(queueLock);
            planQueue.push({plan.first, plan.second});
        }
        currAngles = plan.first->angles;
        ballStart = plan.second;
    }
}

/*
* execution thread reads from the current plan and updates the world state,
* adds the slowdown at the end of the plan
*/
void executionThread() {
    Node *node;
    pair<Node*, BallState> out;
    AngleState angles;
    double t0, tf;
    stack<AngleState> angleStack;
    bool wait, updateBall = false;
    while (true) {
        wait = false;
        if (angleStack.empty()) {
            // finished this path, pop the next from the queue
            {
                lock_guard<mutex> lock(queueLock);
                if (planQueue.empty()) {
                    wait = true;
                } else {
                    out = planQueue.front();
                    planQueue.pop();
                    node = out.first;
                    while (node->bp != nullptr) {
                        angleStack.push(node->angles);
                        node = node->bp;
                    }
                    angleStack.push(node->angles);
                    updateBall = true;
                }
            }
        }
        if (!wait) {
            angles = angleStack.top();
            angleStack.pop();
        }
        {
            lock_guard<mutex> lock(worldLock);
            tf = world.angleEnd.t;
            t0 = world.time;
        }
        this_thread::sleep_for(chrono::milliseconds(int((tf - t0)*1000)));
        {
            lock_guard<mutex> lock(worldLock);
            world.angleStart = world.angleEnd;
            if (!wait) {
                world.angleEnd = angles;
            }
            world.ballStart = world.ballNext;
            if (updateBall && angleStack.empty()) {
                cout << (out.second.t == angles.t) << endl;
                world.ballNext = out.second;
                updateBall = false;
            }
        }
    }
}

// returns rectangle to display for this link from p1 to p2 with thickness
sf::RectangleShape makeLink(sf::Vector2f p1, sf::Vector2f p2, float thickness) {
    sf::Vector2f diff = p2 - p1;
    float length = sqrt(diff.x * diff.x + diff.y * diff.y);
    float angle = atan2(diff.y, diff.x) * 180.0f / PI;

    sf::RectangleShape rect(sf::Vector2f(length, thickness));
    rect.setFillColor(sf::Color::White);

    rect.setOrigin(0, thickness / 2.0f); // center vertically
    rect.setPosition(p1);
    rect.setRotation(angle);

    return rect;
}

/*
* visualizer thread runs the graphics and reads from world state 
* periodically to update graphics
*/
void visualizerThread() {
    double tnow, effectorSize;
    WorldState snapshot;
    sf::RenderWindow window(sf::VideoMode(800, 600), "Arm Visualizer");
    window.setFramerateLimit(60);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        tnow = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()).count() 
            / 1000.0 - startTime;
        {
            lock_guard<mutex> lock(worldLock);
            world.time = tnow;
            snapshot = world;
        }
        BallState b = ball(snapshot.ballStart, snapshot.time);

        EffectorState effector = anglesToEffector(arm(snapshot.angleStart, snapshot.angleEnd, snapshot.time));

        sf::Vector2f base(VISX(BASEX), VISY(BASEY));
        sf::Vector2f joint(VISX(effector.xmid), VISY(effector.ymid));
        sf::Vector2f end(VISX(effector.x), VISY(effector.y));

        window.clear(sf::Color::Black);

        sf::RectangleShape baseRect(sf::Vector2f(50.0, 4.0f));
        baseRect.setFillColor(sf::Color::Blue);
        baseRect.setOrigin(25.0f, 2.0f);
        baseRect.setPosition(base);

        auto link1 = makeLink(base, joint, 4.0f);
        auto link2 = makeLink(joint, end, 4.0f);
        
        window.draw(link1);
        window.draw(link2);
        window.draw(baseRect);

        effectorSize = EFFECTORRADIUS * WINDOWSCALE;
        sf::CircleShape efffector(effectorSize);
        efffector.setFillColor(sf::Color::Green);
        efffector.setOrigin(effectorSize, effectorSize);
        efffector.setPosition(end.x, end.y);
        window.draw(efffector);

        sf::CircleShape circle(BALLRADIUS*WINDOWSCALE);
        circle.setFillColor(sf::Color::Red);
        circle.setOrigin(BALLRADIUS*WINDOWSCALE, BALLRADIUS*WINDOWSCALE);
        circle.setPosition(VISX(b.x), VISY(b.y));
        window.draw(circle);

        window.display();
    }
}

// run planner, execution, and visualizer simultaneously
int main(int argc, char** argv) {
    cout << "Initializing...\n";
    world.angleStart = {PI/4, PI/2, 0, 0, 0};
    world.angleEnd = {PI/4, PI/2, 0, 0, 2};
    world.effector = anglesToEffector(world.angleStart);
    world.ballStart = {world.effector.x, world.effector.y, 0, 10, 2};
    world.ballNext = world.ballStart;
    world.time = 0.0;
    thread execution(executionThread);
    thread planner(plannerThread);
    visualizerThread(); // run in main thread

    execution.join();
    planner.join();

    return 0;
}