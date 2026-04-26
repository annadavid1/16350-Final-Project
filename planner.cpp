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
#include <algorithm>
#include <fstream>


using namespace std;
using namespace sf;

#define G 9.8
#define PI M_PI
#define DOUBLEMAX numeric_limits<double>::max()

#define WINDOWX 1100
#define WINDOWY 900
#define WINDOWSCALE 40
#define FRAMERATE 60
#define FONTSIZE 18

#define VIS(x) WINDOWSCALE * x
#define VISX(x) VIS(x) + WINDOWX/2
#define VISY(y) WINDOWY - VIS(y) - 150
#define SEC(x) x / 1000.0

#define BALLRADIUS 0.35

#define BASEX1 -7.5
#define BASEY1 0.0
#define BASEX2 7.5
#define BASEY2 0.0
#define LINK1 3.0
#define LINK2 3.0

#define BASEWIDTH 0.75
#define BASEHEIGHT 0.08
#define LINKWIDTH 0.08
#define EFFECTORWIDTH 0.5
#define EFFECTORHEIGHT 0.08

#define TARGETX1 BASEX2
#define TARGETY1 BASEY2
#define TARGETX2 BASEX1
#define TARGETY2 BASEY1
#define TARGETRADIUS (LINK1 + 3.0 * LINK2 / 4.0)
#define MINHEIGHT BASEY1 + LINK1 + LINK2 * (4.5 + max(0, NUMBALLS-4) * 0.75)
#define MAXHEIGHT BASEY1 + LINK1 + LINK2 * (5.5 + max(0, NUMBALLS-4) * 0.75)

#define MAXOMEGA1 (2.0 + max(0, NUMBALLS-4) * 0.3)*PI
#define MAXOMEGA2 (2.0 + max(0, NUMBALLS-4) * 0.3)*PI
#define MAXACC1 (3.5 + max(0, NUMBALLS-4) * 0.35)*PI
#define MAXACC2 (3.5 + max(0, NUMBALLS-4) * 0.35)*PI
#define TIMEWEIGHT 1.0

#define TMAX 1
#define DT 0.1
#define DA MAXACC1 / 16

// epsilon for extend in RRT
#define EPSILON PI/16.0

#define SAMPLETIME 2.2
#define FIRSTTHROW 2.5*NUMBALLS + 5.0*max(0,NUMBALLS-4)

// how far ahead we can plan
#define PLANAHEAD 50 + NUMBALLS*50
#define MAXELAPSE 1.0 + 0.5*(NUMBALLS > 4)
#define ELAPSEMULT 1.25

// extend outputs
#define TRAPPED 0
#define ADVANCED 1
#define REACHED 2


// initialize random distribution
random_device rd;
default_random_engine gen(rd());
uniform_real_distribution<> rand_omega1(-MAXOMEGA1, MAXOMEGA1);
uniform_real_distribution<> rand_omega2(-MAXOMEGA2, MAXOMEGA2);


struct BallState
{
    double x;
    double y;
    double vx;
    double vy;
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
    double basex;
    double basey;
    double targetx;
    double targety;
    double theta1;
    double theta2;
    double omega1;
    double omega2;
    double t;
};

struct EffectorState
{
    double basex;
    double basey;
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
    AngleState angleStart[2];
    AngleState angleEnd[2];
    EffectorState effector;
    BallState* ballStart;
    BallState* ballNext;
    int executed;
    int executedNext;
    int planned;
    float time;
};

struct Node
{
    AngleState angles;
    Node *bp;
    Node *left; // left subtree in kd-tree
    Node *right; // right subtree in kd-tree

    Node(AngleState angles) : angles(angles), bp(nullptr), left(nullptr), right(nullptr) {}
    Node(AngleState angles, Node *bp) : angles(angles), bp(bp), left(nullptr), right(nullptr) {}
};

struct RRTOut
{
    Node* node;
    BallState newBall;
    bool done;
};

// tree storing root and list of nodes
struct Tree {
	Node *root; // start pos

	// tree initialized with just root
	Tree(Node *root) : root(root) {};
};

WorldState world;
mutex worldLock;

deque<Node*> planQueue1;
deque<Node*> planQueue2;

deque<BallState>* ballQueue;
mutex queueLock;

Color ballColors[5] = {Color::Yellow, Color::Magenta, Color::Cyan, Color::Red, 
                       Color::Green};

vector<double> times;
vector<double> xs;
vector<double> ys;


double startTime;

std::ofstream resfile;


// get position of ball at time t that started at start
static BallState ball(BallState start, double t) {
    double dt = t - start.t;
    double vx = start.vx;
    double vy = start.vy - G * dt;
    double x = start.x + vx * dt;
    double y = start.y + start.vy * dt - G * dt * dt / 2;
    return {x, y, vx, vy, t};
}

// get position of arm at time t moving from start to end
static AngleState arm(AngleState start, AngleState end, double t) {
    double tf = end.t - start.t;
    double dt = t - start.t;
    double a1 = (end.omega1 - start.omega1)/tf;
    double omega1 = start.omega1 + a1*dt;
    double theta1 = start.theta1 + start.omega1*dt + a1*dt*dt*0.5;
    double a2 = (end.omega2 - start.omega2)/tf;
    double omega2 = start.omega2 + a2*dt;
    double theta2 = start.theta2 + start.omega2*dt + a2*dt*dt*0.5;
    return {start.basex, start.basey, start.targetx, start.targety, 
            theta1, theta2, omega1, omega2, t};
}

// get position of arm at time t moving from start to end
static AngleState arm(AngleState start, OmegaState end, double t) {
    double tf = end.t - start.t;
    double dt = t - start.t;
    double a1 = (end.omega1 - start.omega1)/tf;
    double omega1 = start.omega1 + a1*dt;
    double theta1 = start.theta1 + start.omega1*dt + a1*dt*dt*0.5;
    double a2 = (end.omega2 - start.omega2)/tf;
    double omega2 = start.omega2 + a2*dt;
    double theta2 = start.theta2 + start.omega2*dt + a2*dt*dt*0.5;
    return {start.basex, start.basey, start.targetx, start.targety, theta1, 
            theta2, omega1, omega2, t};
}

// get end effector state from angle state
static EffectorState anglesToEffector(AngleState angles) {
    double theta12 = angles.theta1+angles.theta2;
    double omega12 = angles.omega1 + angles.omega2;
    double cos1 = cos(angles.theta1);
    double sin1 = sin(angles.theta1);
    double cos2 = cos(theta12);
    double sin2 = sin(theta12);
    double xmid = angles.basex + LINK1 * cos1;
    double x = xmid + LINK2 * cos2;
    double ymid = angles.basey + LINK1 * sin1;
    double y = ymid + LINK2 * sin2;
    double vx = -LINK1 * sin1 * angles.omega1 - LINK2 
        * sin2 * omega12;
    double vy = LINK1 * cos1 * angles.omega1 + LINK2 
        * cos2 * omega12;
    return {angles.basex, angles.basey, xmid, ymid, x, y, vx, vy, angles.t};
}

// checks if effector will catch ball and throw will put ball back within target
static pair<bool, BallState> isGoalConfig(
    AngleState angles, 
    BallState currBall
) {
    EffectorState effector = anglesToEffector(angles);
    if (effector.vy <= 0) {
        return {false, {0,0,0,0,0}};
    }
    double t, t0, xdist, ydist, dist, x0, y0, vx, vy, a, b, c, d, e, A, B, C, D,
        p, q, u, v, r, phi, m, minT, maxDist, res;
    xdist = currBall.x - effector.x;
    ydist = currBall.y - BALLRADIUS - effector.y;
    if (abs(xdist) > EFFECTORWIDTH*0.5 || abs(ydist) > EFFECTORHEIGHT*0.5) {
        return {false, {0,0,0,0,0}};
    }
    int i;
    double ts[3];
    x0 = currBall.x - angles.targetx;
    y0 = currBall.y - angles.targety;
    vx = effector.vx;
    vy = effector.vy;
    t0 = effector.t;

    if (y0 + vy*vy/(2*G) < MINHEIGHT - angles.targety || y0 + vy*vy/(2*G) > MAXHEIGHT - angles.targety) {
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
        u = cbrt(-q*0.5 + sqrt(D));
        v = cbrt(-q*0.5 - sqrt(D));
        t = u + v - A/3.0;
        if (t < minT) {
            return {false, {0,0,0,0,0}};
        }
        xdist = x0+vx*t;
        ydist = y0+vy*t-G*t*t*0.5;
        if (xdist*xdist + ydist*ydist <= TARGETRADIUS*TARGETRADIUS) {
            return {true, 
                    {x0 + angles.targetx, y0 + angles.targety, vx, vy, t0}};
        }
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
        res = a*t*t*t*t + b*t*t*t + c*t*t + d*t + e;
        if (res <= 0) {
            return {true, 
                    {x0 + angles.targetx, y0 + angles.targety, vx, vy, t0}};
        }
    }
    return {false, {0,0,0,0,0}};
}

static bool solveCatchForTime(
    const AngleState& start,
    const BallState& ballStart,
    double dt,
    AngleState& outAngles,
    BallState& outBall
) {
    BallState ballt = ball(ballStart, start.t + dt);
    double x, y;

    // desired catch point
    double midx = ballt.x;
    double midy = ballt.y - BALLRADIUS;

    for (double xchange : xs) {
        x = midx + xchange;
        for (double ychange : ys) {
            y = midy + ychange;
            double dx = x - start.basex;
            double dy = y - start.basey;
            double r2 = dx*dx + dy*dy;

            double c2 = (r2 - LINK1*LINK1 - LINK2*LINK2) / (2.0*LINK1*LINK2);
            if (c2 < -1.0 || c2 > 1.0) continue;

            // two IK branches
            double s2a = sqrt(max(0.0, 1.0 - c2*c2));
            double s2b = -s2a;

            for (double s2 : {s2a, s2b}) {
                double theta2p = atan2(s2, c2);
                while (theta2p - start.theta2 > PI) theta2p -= 2.0 * PI;
                while (start.theta2 - theta2p > PI) theta2p += 2.0 * PI;
                double k1 = LINK1 + LINK2 * c2;
                double k2 = LINK2 * s2;
                double theta1p = atan2(dy, dx) - atan2(k2, k1);
                while (theta1p - start.theta1 > PI) theta1p -= 2.0 * PI;
                while (start.theta1 - theta1p > PI) theta1p += 2.0 * PI;

                double a1 = 2.0 * (theta1p - start.theta1 - start.omega1*dt) / (dt*dt);
                double a2 = 2.0 * (theta2p - start.theta2 - start.omega2*dt) / (dt*dt);

                if (a1 < -MAXACC1 || a1 > MAXACC1 || a2 < -MAXACC2 || a2 > MAXACC2) {
                    continue;
                }

                double omega1p = start.omega1 + a1*dt;
                double omega2p = start.omega2 + a2*dt;

                if (omega1p < -MAXOMEGA1 || omega1p > MAXOMEGA1 ||
                    omega2p < -MAXOMEGA2 || omega2p > MAXOMEGA2) {
                    continue;
                }

                AngleState candidate = {
                    start.basex, start.basey, start.targetx, start.targety,
                    theta1p, theta2p, omega1p, omega2p, start.t + dt
                };

                auto goal = isGoalConfig(candidate, ballt);
                if (goal.first) {
                    outAngles = candidate;
                    outBall = goal.second;
                    return true;
                }
            }
        }
    }

    return false;
}

// checks if a state can be created from this state which is a goal state
static pair<bool, BallState> isGoalConfigFrom(
    AngleState angles, 
    BallState ballStart, 
    AngleState& anglesp
) {
    BallState newBall;
    for (double t : times) {
        if (solveCatchForTime(angles, ballStart, t, anglesp, newBall)) {
            return {true, newBall};
        }
    }
    return {false, {0,0,0,0,0}};
}

// Euclidean distance between angle and omega configs
static double dist(AngleState start, OmegaState end) {
	double omega1 = start.omega1 - end.omega1;
    double omega2 = start.omega2 - end.omega2;
    double t = (start.t - end.t)*TIMEWEIGHT; // scale time down so not overpower
	return sqrt(omega1*omega1 + omega2*omega2 + t*t);
}

// check end after start and acceleration less than max
static bool canConnect(AngleState start, OmegaState end);

// recursively inserts node into kd-tree
static void insert(Node *t, Node *v, int depth) {
	int k = depth % 3;
    double tval, vval;
    if (k == 0) {
        tval = t->angles.omega1;
        vval = v->angles.omega1;
    } else if (k == 1) {
        tval = t->angles.omega2;
        vval = v->angles.omega2;
    } else {
        tval = t->angles.t * TIMEWEIGHT;
        vval = v->angles.t * TIMEWEIGHT;
    }
	if (vval <= tval) {
		// less than current, move left
		if (t->left == nullptr) {
			// found leaf
			t->left = v;
		} else {
			insert(t->left, v, depth+1);
		}
	} else {
		// greater than current, move right
		if (t->right == nullptr) {
			// found leaf
			t->right = v;
		} else {
			insert(t->right, v, depth+1);
		}
	}
}

// finds nearest neighbor using kd-tree
static Node* nearest(Node* t, OmegaState v, int depth) {
    if (t == nullptr) {
        return nullptr;
    }
	int k = depth % 3;
    double tval, vval;
    if (k == 0) {
        tval = t->angles.omega1;
        vval = v.omega1;
    } else if (k == 1) {
        tval = t->angles.omega2;
        vval = v.omega2;
    } else {
        tval = t->angles.t * TIMEWEIGHT;
        vval = v.t * TIMEWEIGHT;
    }
	double d, mindist, dt, bestDist = DOUBLEMAX;
	Node *best = nullptr, *n, *n2, *nearChild, *farChild;

	dt = dist(t->angles, v);
	mindist = fabs(vval-tval); // min dist from any node on far side

    if (canConnect(t->angles, v)) {
        best = t;
        bestDist = dt;
    }

    if (vval <= tval) {
        nearChild = t->left;
        farChild = t->right;
    } else {
        nearChild = t->right;
        farChild = t->left;
    }

    if (nearChild != nullptr) {
        n = nearest(nearChild, v, depth+1);
        if (n != nullptr) {
            d = dist(n->angles, v);
            if (d < bestDist) {
                best = n;
                bestDist = d;
            }
        }
    }

    if (farChild != nullptr && bestDist > mindist) {
        n2 = nearest(farChild, v, depth+1);
        if (n2 != nullptr) {
            d = dist(n2->angles, v);
            if (d < bestDist) {
                best = n2;
            }
        }
    }

	return best;
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

// for RRT, sample random state and extend tree towards it, return new state
static pair<int, Node*> extend(Tree* t, BallState ballStart, double startTime) {
    // sample from startTime to SAMPLETIME * (time to max ball height)
    OmegaState qrand = sample(startTime, 
        ballStart.t + ballStart.vy / G * SAMPLETIME);
    Node *q, *qnew, *qmin = nullptr;
    double d, scale, minDist = DOUBLEMAX;
    OmegaState qext;
    qmin = nearest(t->root, qrand, 0);
    if (qmin == nullptr) {
        return {TRAPPED, nullptr};
    }
    minDist = dist(qmin->angles, qrand);
    if (minDist == 0) {
        return {TRAPPED, nullptr};
    }
    scale = EPSILON / minDist;
    if (scale >= 1) {
        qnew = new Node(arm(qmin->angles, qrand, qrand.t), qmin);
        insert(t->root, qnew, 0);
        return {REACHED, qnew};
    }

    qext = {
        qmin->angles.omega1 + (qrand.omega1 - qmin->angles.omega1) * scale,
        qmin->angles.omega2 + (qrand.omega2 - qmin->angles.omega2) * scale,
        qmin->angles.t + (qrand.t - qmin->angles.t) * scale
    };
    qnew = new Node(arm(qmin->angles, qext, qext.t), qmin);
    insert(t->root, qnew, 0);
    return {ADVANCED, qnew};
}

// RRT algorithm randomly sampling omegas and time to extend
static RRTOut runRRT(
    AngleState startAngle,
    BallState ballStart,
    double elapseMult
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

    double rrtStart = SEC(chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch()).count());

	// main loop, goes until any path found
	while (1) {
		// get random config and extend t towards it
        if (SEC(chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch()).count()) > rrtStart+MAXELAPSE*elapseMult) {
            return {{}, {}, false};
        }
		ext = extend(t, ballStart, startAngle.t);
		s = ext.first; qnew = ext.second;
		if (s != TRAPPED) {
            // checks if can create a goal config after this state
            isGoal = isGoalConfigFrom(qnew->angles, ballStart, anglesp);
            if (isGoal.first) {
                // reached goal
                qnew = new Node(anglesp, qnew);
                // no need to insert into kd-tree since done with tree now
                return {qnew, isGoal.second, true};
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
    RRTOut plan1, plan2;
    WorldState snapshot;
    int planned, executed;
    double time;
    {
        lock_guard<mutex> lock(worldLock);
        snapshot = world;
    }
    AngleState currAngles1 = snapshot.angleEnd[0];
    AngleState currAngles2 = snapshot.angleEnd[1];
    BallState* ballWorld = snapshot.ballStart;
    BallState ballStart[NUMBALLS];
    for (int i = 0; i < NUMBALLS; i++) {
        ballStart[i] = ballWorld[i];
    }
    int ballIndex = 0;
    int backtrack = 1;
    int backtrackPlan = -1;
    double elapseMult = 1;
    while (1)
    {
        {
            lock_guard<mutex> lock(worldLock);
            snapshot = world;
        }
        if (snapshot.planned >= snapshot.executed + PLANAHEAD) {
            this_thread::sleep_for(chrono::milliseconds(1000));
            continue;
        }

        // left arm
        while (true) {
            {
                lock_guard<mutex> lock(queueLock);
                if (!planQueue1.empty()) {
                    currAngles1 = planQueue1.back()->angles;
                }
                if (!ballQueue[ballIndex].empty()) {
                    ballStart[ballIndex] = ballQueue[ballIndex].back();
                }
            }
            plan1 = runRRT(currAngles1, ballStart[ballIndex], elapseMult);
            if (plan1.done) {
                break;
            } else {
                shuffle(times.begin(), times.end(), gen);
                shuffle(xs.begin(), xs.end(), gen);
                shuffle(ys.begin(), ys.end(), gen);
                {
                    lock_guard<mutex> lock(worldLock);
                    if (backtrackPlan == -1) {
                        backtrackPlan = world.planned;
                    }
                    if (world.planned - backtrack*2 < world.executed + 2) {
                        backtrack = (double)(world.planned - world.executed - 2)*0.5;
                    }
                    for (int i = 0; i < backtrack; i++) {
                        ballIndex = ballIndex > 0 ? ballIndex-1 : NUMBALLS-1;
                        {
                            lock_guard<mutex> lock(queueLock);
                            planQueue1.pop_back();
                            planQueue2.pop_back();
                            ballQueue[ballIndex].pop_back();
                            ballQueue[ballIndex].pop_back();
                        }
                    }
                }
                {
                    lock_guard<mutex> lock(worldLock);
                    world.planned -= backtrack * 2;
                    planned = world.planned;
                    time = world.time;
                    executed = world.executedNext;
                }
                resfile << time << ", " << planned << ", " << executed << endl;
                backtrack++;
                elapseMult *= ELAPSEMULT;
            }
        }
        {
            lock_guard<mutex> lock(queueLock);
            planQueue1.push_back(plan1.node);
            ballQueue[ballIndex].push_back(plan1.newBall);
        }
        {
            lock_guard<mutex> lock(worldLock);
            world.planned += 1;
            if (backtrackPlan < world.planned) {
                backtrackPlan = -1;
                backtrack = 1;
                elapseMult = 1;
            }
            planned = world.planned;
            time = world.time;
            executed = world.executedNext;
        }
        resfile << time << ", " << planned << ", " << executed << endl;
        currAngles1 = plan1.node->angles;
        ballStart[ballIndex] = plan1.newBall;

        // right arm
        while (true)
        {
            {
                lock_guard<mutex> lock(queueLock);
                if (!planQueue2.empty()) {
                    currAngles2 = planQueue2.back()->angles;
                }
                if (!ballQueue[ballIndex].empty()) {
                    ballStart[ballIndex] = ballQueue[ballIndex].back();
                }
            }
            plan2 = runRRT(currAngles2, ballStart[ballIndex], elapseMult);
            if (plan2.done) {
                break;
            } else {
                shuffle(times.begin(), times.end(), gen);
                shuffle(xs.begin(), xs.end(), gen);
                shuffle(ys.begin(), ys.end(), gen);
                {
                    lock_guard<mutex> lock(worldLock);
                    if (backtrackPlan == -1) {
                        backtrackPlan = world.planned;
                    }
                    if (world.planned - backtrack*2 < world.executed + 2) {
                        backtrack = (double)(world.planned - world.executed - 2)*0.5;
                    }
                    for (int i = 0; i < backtrack; i++) {
                        {
                            lock_guard<mutex> lock(queueLock);
                            planQueue1.pop_back();
                            planQueue2.pop_back();
                            ballQueue[ballIndex].pop_back();
                            ballIndex = ballIndex > 0 ? ballIndex-1 : NUMBALLS-1;
                            ballQueue[ballIndex].pop_back();
                        }
                    }
                }
                {
                    lock_guard<mutex> lock(worldLock);
                    world.planned -= backtrack * 2;
                    planned = world.planned;
                    time = world.time;
                    executed = world.executedNext;
                }
                resfile << time << ", " << planned << ", " << executed << endl;
                backtrack++;
                elapseMult *= ELAPSEMULT;
            }
        }
        
        {
            lock_guard<mutex> lock(queueLock);
            planQueue2.push_back(plan2.node);
            ballQueue[ballIndex].push_back(plan2.newBall);
        }
        {
            lock_guard<mutex> lock(worldLock);
            world.planned += 1;
            if (backtrackPlan < world.planned) {
                backtrackPlan = -1;
                backtrack = 1;
                elapseMult = 1;
            }
            planned = world.planned;
            time = world.time;
            executed = world.executedNext;
        }
        resfile << time << ", " << planned << ", " << executed << endl;
        currAngles2 = plan2.node->angles;
        ballStart[ballIndex] = plan2.newBall;
        ballIndex++;
        if (ballIndex >= NUMBALLS) {
            ballIndex = 0;
        }
    }
}

/*
* execution thread reads from the current plan and updates the world state,
* adds the slowdown at the end of the plan
*/
void executionThread() {
    Node *node, *prev;
    AngleState angles1, angles2;
    double t0, tf1, tf2;
    stack<AngleState> angleStack1, angleStack2;
    bool wait1, wait2, updateBall1 = false, updateBall2 = false, 
        nextBall1 = false, nextBall2 = false;
    int ballIndex1 = 0, ballIndex2 = 0;
    while (true) {
        wait1 = false;
        wait2 = false;
        {
            lock_guard<mutex> lock(worldLock);
            tf1 = world.angleEnd[0].t;
            tf2 = world.angleEnd[1].t;
            t0 = world.time;
        }

        if (tf1 < tf2) {
            // left arm
            if (angleStack1.empty()) {
                // finished this path, pop the next from the queue
                {
                    lock_guard<mutex> lock(queueLock);
                    if (planQueue1.empty()) {
                        if (t0 >= FIRSTTHROW) {
                            resfile.close();
                            cout << "No plans left\n";
                            exit(0);
                        }
                        wait1 = true;
                    } else {
                        node = planQueue1.front();
                        planQueue1.pop_front();
                        while (node->bp != nullptr) {
                            angleStack1.push(node->angles);
                            prev = node;
                            node = node->bp;
                            delete prev;
                        }
                        angleStack1.push(node->angles);
                        delete node;
                        updateBall1 = true;
                    }
                }
            }
            if (!wait1) {
                angles1 = angleStack1.top();
                angleStack1.pop();
            }
            {
                lock_guard<mutex> lock(worldLock);
                tf1 = world.angleEnd[0].t;
                t0 = world.time;
            }
            this_thread::sleep_for(chrono::milliseconds(int((tf1 - t0)*1000)));
            {
                lock_guard<mutex> lock(worldLock);
                world.angleStart[0] = world.angleEnd[0];
                if (!wait1) {
                    world.angleEnd[0] = angles1;
                }
                if (nextBall1) {
                    world.ballStart[ballIndex1] = world.ballNext[ballIndex1];
                    world.executed = world.executedNext;
                    nextBall1 = false;
                    ballIndex1++;
                    if (ballIndex1 >= NUMBALLS) {
                        ballIndex1 = 0;
                    }
                }
                lock_guard<mutex> lock1(queueLock);
                if (updateBall1 && angleStack1.empty()) {
                    if (ballQueue[ballIndex1].empty()) {
                        cout << "No plans left\n";
                        exit(0);
                    }
                    world.executedNext += 1;
                    world.ballNext[ballIndex1] = ballQueue[ballIndex1].front();
                    ballQueue[ballIndex1].pop_front();
                    updateBall1 = false;
                    nextBall1 = true;
                    resfile << world.time << ", " << world.planned << ", " << world.executedNext << endl;
                }
            }
        } else {
            // right arm
            if (angleStack2.empty()) {
                // finished this path, pop the next from the queue
                {
                    lock_guard<mutex> lock(queueLock);
                    if (planQueue2.empty()) {
                        wait2 = true;
                        if (t0 >= 10) {
                            resfile.close();
                            cout << "No plans left\n";
                            exit(0);
                        }
                    } else {
                        node = planQueue2.front();
                        planQueue2.pop_front();
                        while (node->bp != nullptr) {
                            angleStack2.push(node->angles);
                            prev = node;
                            node = node->bp;
                            delete prev;
                        }
                        angleStack2.push(node->angles);
                        delete node;
                        updateBall2 = true;
                    }
                }
            }
            if (!wait2) {
                angles2 = angleStack2.top();
                angleStack2.pop();
            }
            {
                lock_guard<mutex> lock(worldLock);
                tf2 = world.angleEnd[1].t;
                t0 = world.time;
            }
            this_thread::sleep_for(chrono::milliseconds(int((tf2 - t0)*1000)));
            {
                lock_guard<mutex> lock(worldLock);
                world.angleStart[1] = world.angleEnd[1];
                if (!wait2) {
                    world.angleEnd[1] = angles2;
                }
                if (nextBall2) {
                    world.ballStart[ballIndex2] = world.ballNext[ballIndex2];
                    world.executed = world.executedNext;
                    nextBall2 = false;
                    ballIndex2++;
                    if (ballIndex2 >= NUMBALLS) {
                        ballIndex2 = 0;
                    }
                }
                lock_guard<mutex> lock1(queueLock);
                if (updateBall2 && angleStack2.empty()) {
                    if (ballQueue[ballIndex2].empty()) {
                        cout << "No plans left\n";
                        exit(0);
                    }
                    world.executedNext += 1;
                    world.ballNext[ballIndex2] = ballQueue[ballIndex2].front();
                    ballQueue[ballIndex2].pop_front();
                    updateBall2 = false;
                    nextBall2 = true;
                    resfile << world.time << ", " << world.planned << ", " << world.executedNext << endl;
                }
            }
        }
    }
}

// returns rectangle to display for this link from p1 to p2 with thickness
RectangleShape makeLink(Vector2f p1, Vector2f p2, float thickness) {
    Vector2f diff = p2 - p1;
    float length = sqrt(diff.x * diff.x + diff.y * diff.y);
    float angle = atan2(diff.y, diff.x) * 180.0f / PI;

    RectangleShape rect(Vector2f(length, thickness));
    rect.setFillColor(Color::White);

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
    float padding = 10.0f;
    WorldState snapshot;
    Font font;
    Text executedText, plannedText, timeText;
    Event event;
    BallState b;
    EffectorState eff1, eff2;
    FloatRect bounds1, bounds2, bounds3;
    
    if (!font.loadFromFile("Arial.ttf")) {
        cout << "Could not load font\n";
        exit(0);
    }

    RenderWindow window(VideoMode(WINDOWX, WINDOWY), "Arm Visualizer");
    window.setFramerateLimit(FRAMERATE);

    executedText.setFont(font);
    executedText.setCharacterSize(FONTSIZE);
    executedText.setFillColor(Color::White);

    plannedText.setFont(font);
    plannedText.setCharacterSize(FONTSIZE);
    plannedText.setFillColor(Color::White);

    timeText.setFont(font);
    timeText.setCharacterSize(FONTSIZE);
    timeText.setFillColor(Color::White);

    while (window.isOpen()) {
        while (window.pollEvent(event)) {
            if (event.type == Event::Closed) {
                window.close();
            }
        }

        tnow = SEC(chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch()).count())-startTime;
        {
            lock_guard<mutex> lock(worldLock);
            world.time = tnow;
            snapshot = world;
        }

        // left arm
        eff1 = anglesToEffector(
            arm(snapshot.angleStart[0], snapshot.angleEnd[0], snapshot.time));

        Vector2f base1(VISX(BASEX1), VISY(BASEY1));
        Vector2f joint1(VISX(eff1.xmid), VISY(eff1.ymid));
        Vector2f end1(VISX(eff1.x), VISY(eff1.y));

        window.clear(Color::Black);

        RectangleShape baseRect1(Vector2f(VIS(BASEWIDTH), VIS(BASEHEIGHT)));
        baseRect1.setFillColor(Color::Blue);
        baseRect1.setOrigin(VIS(BASEWIDTH*0.5f), VIS(BASEHEIGHT*0.5f));
        baseRect1.setPosition(base1);

        auto link11 = makeLink(base1, joint1, VIS(LINKWIDTH));
        auto link21 = makeLink(joint1, end1, VIS(LINKWIDTH));
        
        window.draw(link11);
        window.draw(link21);
        window.draw(baseRect1);

        RectangleShape effector1(Vector2f(VIS(EFFECTORWIDTH), 
            VIS(EFFECTORHEIGHT)));
        effector1.setFillColor(Color::Green);
        effector1.setOrigin(VIS(EFFECTORWIDTH*0.5f), VIS(EFFECTORHEIGHT*0.5f));
        effector1.setPosition(end1);

        window.draw(effector1);

        // right arm
        eff2 = anglesToEffector(
            arm(snapshot.angleStart[1], snapshot.angleEnd[1], snapshot.time));

        Vector2f base2(VISX(BASEX2), VISY(BASEY2));
        Vector2f joint2(VISX(eff2.xmid), VISY(eff2.ymid));
        Vector2f end2(VISX(eff2.x), VISY(eff2.y));

        RectangleShape baseRect2(Vector2f(VIS(BASEWIDTH), VIS(BASEHEIGHT)));
        baseRect2.setFillColor(Color::Blue);
        baseRect2.setOrigin(VIS(BASEWIDTH*0.5f), VIS(BASEHEIGHT*0.5f));
        baseRect2.setPosition(base2);

        auto link12 = makeLink(base2, joint2, VIS(LINKWIDTH));
        auto link22 = makeLink(joint2, end2, VIS(LINKWIDTH));
        
        window.draw(link12);
        window.draw(link22);
        window.draw(baseRect2);

        RectangleShape effector2(Vector2f(VIS(EFFECTORWIDTH), 
            VIS(EFFECTORHEIGHT)));
        effector2.setFillColor(Color::Green);
        effector2.setOrigin(VIS(EFFECTORWIDTH*0.5f), VIS(EFFECTORHEIGHT*0.5f));
        effector2.setPosition(end2);

        window.draw(effector2);

        // ball
        for (int i = 0; i < NUMBALLS; i++) {
            b = ball(snapshot.ballStart[i], snapshot.time);
            CircleShape circle(VIS(BALLRADIUS));
            circle.setFillColor(ballColors[i%5]);
            circle.setOrigin(VIS(BALLRADIUS), VIS(BALLRADIUS));
            circle.setPosition(VISX(b.x), VISY(b.y));
            window.draw(circle);
        }

        // text
        executedText.setString("Executed Catches: " 
            + to_string(snapshot.executed));
        plannedText.setString("Planned Catches: " 
            + to_string(snapshot.planned));
        timeText.setString("Time (s): " 
            + to_string((int)(snapshot.time)));

        bounds1 = executedText.getLocalBounds();
        executedText.setPosition(
            window.getSize().x - bounds1.width - padding,
            padding
        );

        bounds2 = plannedText.getLocalBounds();
        plannedText.setPosition(
            window.getSize().x - bounds2.width - padding,
            padding + bounds1.height + 10
        );

        bounds3 = timeText.getLocalBounds();
        timeText.setPosition(
            window.getSize().x - bounds3.width - padding,
            padding + bounds2.height + bounds1.height + 20
        );

        window.draw(executedText);
        window.draw(plannedText);
        window.draw(timeText);

        window.display();
    }
}

// run planner, execution, and visualizer simultaneously
int main(int argc, char** argv) {
    if (argc != 2) {
        cout << "bad argument\n";
        return 0;
    }
    resfile.open(argv[1]);
    if (resfile.is_open()) {
        resfile << "Time, Plan, Executed\n";
    } else {
        cout << "Couldn't open file\n";
        return 0;
    }
    for (double t = 0; t <= TMAX; t += DT) {
        times.push_back(t);
    }
    for (double x = -EFFECTORWIDTH*0.5; x <= EFFECTORWIDTH*0.5; x += EFFECTORWIDTH/2.0) {
        xs.push_back(x);
    }
    for (double y = -EFFECTORHEIGHT*0.5; y <= EFFECTORHEIGHT*0.5; y += EFFECTORHEIGHT/2.0) {
        ys.push_back(y);
    }
    shuffle(times.begin(), times.end(), gen);
    shuffle(xs.begin(), xs.end(), gen);
    shuffle(ys.begin(), ys.end(), gen);
    ballQueue = new deque<BallState>[NUMBALLS];
    cout << "Initializing...\n";
    world.angleStart[0] = {BASEX1, BASEY1, TARGETX1, TARGETY1, PI*0.25, PI*0.5, 0,0,0};
    world.angleEnd[0] = {BASEX1, BASEY1, TARGETX1, TARGETY1, PI*0.25, PI*0.5, 0, 0, FIRSTTHROW};
    world.angleStart[1] = {BASEX2, BASEY2, TARGETX2, TARGETY2, PI*0.25, PI*0.5, 0,0,0};
    world.angleEnd[1] = {BASEX2, BASEY2, TARGETX2, TARGETY2, PI*0.25, PI*0.5, 0, 0, FIRSTTHROW};
    EffectorState effector = anglesToEffector(world.angleStart[0]);
    world.ballStart = (BallState*)malloc(NUMBALLS * sizeof(BallState));
    world.ballNext = (BallState*)malloc(NUMBALLS * sizeof(BallState));
    for (int i = 0; i < NUMBALLS; i++) {
        world.ballStart[i] = {effector.x, effector.y, 0, 15, FIRSTTHROW+i*8.0/NUMBALLS};
        world.ballNext[i] = world.ballStart[i];
    }
    startTime = SEC(chrono::duration_cast<chrono::milliseconds>(
        chrono::system_clock::now().time_since_epoch()).count());
    world.time = 0.0;
    thread execution(executionThread);
    thread planner(plannerThread);
    cout << "Intial throw in " << FIRSTTHROW << " seconds\n";
    visualizerThread(); // run in main thread
    resfile.close();

    execution.join();
    planner.join();

    return 1;
}
