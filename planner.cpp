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
#define BALLRADIUS 0.15

#define BASEX 0.0
#define BASEY 0.0
#define LINK1 1.0
#define LINK2 1.0

#define EFFECTORRADIUS 0.075

#define TARGETX 0.0
#define TARGETY BASEY + LINK1
#define TARGETRADIUS (LINK2 * 3.0 / 4.0)
#define MINHEIGHT BASEY + LINK1 + 2.0

#define MAXOMEGA1 PI/4.0
#define MAXOMEGA2 PI/4.0
#define MAXACC1 PI/4.0
#define MAXACC2 PI/4.0

// epsilon for extend in RRT
#define EPSILON 2*PI

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

queue<pair<Node*, BallState>> nodeQueue;
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

// check whether a given effector state can successfully catch/throw ball
// static pair<bool, BallState> isGoalConfig(EffectorState effector, BallState ballStart) {
//     BallState b = ball(ballStart, effector.t);
//     int i;
//     double t, t0, xdist, ydist, dist, x0, y0, vx, vy, c, d, D, p, q, u, v, r, phi, m, minT, maxDist;
//     double vys[7], ts[3];
//     xdist = b.x - effector.x;
//     ydist = b.y - effector.y;
//     dist = xdist * xdist + ydist * ydist;
//     maxDist = EFFECTORRADIUS + BALLRADIUS/WINDOWSCALE;
//     if (dist > maxDist * maxDist) {
//         return {false, {0,0,0,0,0}};
//     }
//     x0 = b.x;
//     y0 = b.y;
//     vx = effector.vx;
//     vy = effector.vy;
//     t0 = effector.t;
//     vys[0] = 1;
//     for (i = 1; i < 7; i++) {
//         vys[i] = vys[i-1] * vy;
//     }
//     if (vys[2]/(2*G) < 3.0) {
//         return {false, {0,0,0,0,0}};
//     }
//     minT = vy / G;
//     c = vx*vx + vy*vy - y0*G;
//     d = 2*x0*vx + 2*y0*vy;

//     D = (vys[6]*gs[6] - 2.0/3*vys[4]*gs[6]*c - vys[3]*gs[7]*d 
//         - vy*gs[5]*c*d + vy*gs[7]*c*d + 1.0/4*gs[8]*d*d - 1.0/16*vys[6]*gs[10] 
//         + 1.0/24*vys[4]*gs[10]*c - 4.0/3*vys[2]*gs[6]*c*c - 8.0/27*vys[3]*gs[9])
//         / gs[12];

//     p = (3*vys[2] + 2*c) / gs[2];
//     q = (2*vys[3] + 2*vy*c) / gs[3] + d/gs[2];

//     if (D >= 0) {
//         // one real root
//         u = cbrt(-q/2.0 + sqrt(D));
//         v = cbrt(-q/2.0 - sqrt(D));
//         t = u + v + vy/G;
//         if (t < minT) {
//             return {false, {0,0,0,0,0}};
//         }
//         xdist = x0+vx*t;
//         ydist = y0+vy*t-G*t*t/2;
//         cout << "maybe not goal 3\n";
//         return {xdist*xdist + ydist*ydist <= TARGETRADIUS*TARGETRADIUS, {x0, y0, vx, vy, t0}};
//     }

//     r = sqrt(-p*p*p/27.0);
//     phi = acos(-q/(2.0*r));
//     m = 2 * sqrt(-p/3.0);

//     ts[0] = m * cos(phi/3.0) + vy/G;
//     ts[1] = m * cos((phi + 2.0*PI)/3.0) + vy/G;
//     ts[2] = m * cos((phi + 4.0*PI)/3.0) + vy/G;
//     for (i = 0; i < 3; i++) {
//         t = ts[i];
//         if (t < minT) {
//             continue;
//         }
//         xdist = (x0+vx*t);
//         ydist = (y0+vy*t-G*t*t/2);
//         if (xdist*xdist + ydist*ydist <= TARGETRADIUS*TARGETRADIUS) {
//             cout << "definitely goal 5\n";
//             return {true, {x0, y0, vx, vy, t0}};
//         }
//     }
//     return {false, {0,0,0,0,0}};
// }

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
        cout << "not goal 1\n";
        cout << ydist << "," << xdist << "," << currBall.y << "," << effector.y << "," << effector.t << endl;
        return {false, {0,0,0,0,0}};
    }
    x0 = currBall.x - TARGETX;
    y0 = currBall.y - TARGETY;
    vx = effector.vx;
    vy = effector.vy;
    t0 = effector.t;

    if (y0 + vy*vy/(2*G) < MINHEIGHT - TARGETY) {
        cout << "not goal 2\n";
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
            cout << "not goal 3\n";
            return {false, {0,0,0,0,0}};
        }
        xdist = x0+vx*t;
        ydist = y0+vy*t-G*t*t/2;
        if (xdist*xdist + ydist*ydist <= TARGETRADIUS*TARGETRADIUS) {
            cout << "definitely goal 4\n";
            return {true, {x0 + TARGETX, y0 + TARGETY, vx, vy, t0}};
        }
        cout << "not goal 5\n";
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
            cout << "definitely goal 6\n";
            return {true, {x0 + TARGETX, y0 + TARGETY, vx, vy, t0}};
        }
    }
    cout << "not goal 7\n";
    return {false, {0,0,0,0,0}};
}

// check whether a given effector state can successfully catch/throw ball
static pair<bool, BallState> isGoalConfig(AngleState angles, BallState ballStart) {
    return isGoalConfig(anglesToEffector(angles), ballStart);
}

// Euclidean distance between 2 angle configs
static double dist(AngleState start, AngleState end) {
	double omega1 = start.omega1 - end.omega1;
    double omega2 = start.omega2 - end.omega2;
    double t = start.t - end.t;
	return sqrt(omega1*omega1 + omega2*omega2);
}

// Euclidean distance between angle and omega configs
static double dist(AngleState start, OmegaState end) {
	double omega1 = start.omega1 - end.omega1;
    double omega2 = start.omega2 - end.omega2;
    double t = start.t - end.t;
	return sqrt(omega1*omega1 + omega2*omega2);
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
    return a1 < MAXACC1 && a2 < MAXACC2;
}

// get random angular velocities and time after startT
OmegaState sample(double startT, double endT) {
    double omega1 = rand_omega1(gen);
    double omega2 = rand_omega2(gen);
    uniform_real_distribution<> rand_time(startT, endT);
    double t = rand_time(gen);
    return {omega1, omega2, t};
}

static pair<int, Node*> extend(Tree *t, OmegaState qrand) {
    int n = t->V.size();
    OmegaState qext;
    AngleState qextAngle;
    Node *q, *qnew, *qmin = nullptr;
    double d, scale, minDist = DOUBLEMAX;
    for (int i = 0; i < n; i++) {
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
    scale = EPSILON / minDist;
    if (scale >= 1) {
        qextAngle = arm(qmin->angles, qrand, qrand.t);
        qnew = new Node(qextAngle, qmin);
        return {2, qnew};
    }
    qext = {qmin->angles.omega1 + (qrand.omega1 - qmin->angles.omega1)*scale,
            qmin->angles.omega2 + (qrand.omega2 - qmin->angles.omega2)*scale,
            qmin->angles.t + (qrand.t - qmin->angles.t)*scale};

    qextAngle = arm(qmin->angles, qext, qext.t);
    qnew = new Node(qextAngle, qmin);
    return {1, qnew};
}

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

	// initialize start and goal nodes and trees
	start = new Node(startAngle);
	t = new Tree(start);

	// main loop, goes until any path found
	while (1) {
		// get random config and extend t towards it
		qrand = sample(ballStart.t, ballStart.vy/G*2.25);
		ext = extend(t, qrand);
		s = ext.first; qnew = ext.second;
		if (s != TRAPPED) {
            isGoal = isGoalConfig(qnew->angles, ballStart);
            if (isGoal.first) {
                // reached goal
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
            nodeQueue.push({plan.first, plan.second});
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
    bool wait;
    while (true) {
        wait = false;
        if (angleStack.empty()) {
            // finished this path, pop the next from the queue
            {
                lock_guard<mutex> lock(queueLock);
                if (nodeQueue.empty()) {
                    wait = true;
                } else {
                    out = nodeQueue.front();
                    nodeQueue.pop();
                    node = out.first;
                    while (node->bp != nullptr) {
                        angleStack.push(node->angles);
                        node = node->bp;
                    }
                    angleStack.push(node->angles);
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
            if (angleStack.empty()) {
                world.ballNext = out.second;
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
    world.angleEnd = {PI/4, PI/2, 0, 0, 0};
    world.effector = anglesToEffector(world.angleStart);
    world.ballStart = {world.effector.x, world.effector.y, 0, 10, 0};
    world.ballNext = world.ballStart;
    world.time = startTime;
    thread execution(executionThread);
    thread planner(plannerThread);
    visualizerThread(); // run in main thread

    execution.join();
    planner.join();

    return 0;
}