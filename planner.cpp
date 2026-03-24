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
using namespace sf;

#define G 9.8
#define PI M_PI
#define DOUBLEMAX numeric_limits<double>::max()

#define WINDOWX 800
#define WINDOWY 600
#define WINDOWSCALE 100

#define VISX(x) WINDOWSCALE*x + WINDOWX/2
#define VISY(y) WINDOWY - WINDOWSCALE*y - 100
#define BALLRADIUS 0.2

#define BASEX 0.0
#define BASEY 0.0
#define LINK1 1.5
#define LINK2 1.5

#define EFFECTORRADIUS 0.1

#define TARGETX 0.0
#define TARGETY BASEY + LINK1
#define TARGETRADIUS (LINK2 * 3.0 / 4.0)
#define MINHEIGHT BASEY + LINK1 + LINK2 * 2

#define MAXOMEGA1 PI//2*PI
#define MAXOMEGA2 PI//2*PI
#define MAXACC1 2*PI//4*PI
#define MAXACC2 2*PI//4*PI

#define T_MAX 1
#define DT 0.1
#define DA MAXACC1 / 16

// epsilon for extend in RRT
#define EPSILON PI/8.0

#define SAMPLETIME 2.25

// how far ahead we can plan (up to 5 times as many plans as have been executed)
#define PLANAHEAD 5

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
    int executed;
    int executedNext;
    int planned;
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

double startTime = chrono::duration_cast<chrono::milliseconds>(
    chrono::system_clock::now().time_since_epoch()).count() / 1000.0;


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
    double vx = -LINK1 * sin(angles.theta1) * angles.omega1 - LINK2 
        * sin(angles.theta1 + angles.theta2) * (angles.omega1 + angles.omega2);
    double vy = LINK1 * cos(angles.theta1) * angles.omega1 + LINK2 
        * cos(angles.theta1 + angles.theta2) * (angles.omega1 + angles.omega2);
    return {xmid, ymid, x, y, vx, vy, angles.t};
}

// checks if effector will catch ball and throw will put ball back within target
static pair<bool, BallState> isGoalConfig(
    AngleState angles, 
    BallState ballStart
) {
    EffectorState effector = anglesToEffector(angles);
    if (effector.vy <= 0) {
        return {false, {0,0,0,0,0}};
    }
    BallState currBall = ball(ballStart, effector.t);
    int i;
    double t, t0, xdist, ydist, dist, x0, y0, vx, vy, a, b, c, d, e, A, B, C, D,
        p, q, u, v, r, phi, m, minT, maxDist, res;
    double ts[3];
    xdist = currBall.x - effector.x;
    ydist = currBall.y - effector.y;
    dist = xdist * xdist + ydist * ydist;
    maxDist = EFFECTORRADIUS + BALLRADIUS;
    if (dist > maxDist * maxDist) {
        return {false, {0,0,0,0,0}};
    }
    x0 = currBall.x - TARGETX;
    y0 = currBall.y - TARGETY;
    vx = effector.vx;
    vy = effector.vy;
    t0 = effector.t;

    if (y0 + vy*vy/(2*G) < MINHEIGHT - TARGETY) {
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
            return {false, {0,0,0,0,0}};
        }
        xdist = x0+vx*t;
        ydist = y0+vy*t-G*t*t/2;
        if (xdist*xdist + ydist*ydist <= TARGETRADIUS*TARGETRADIUS) {
            return {true, {x0 + TARGETX, y0 + TARGETY, vx, vy, t0}};
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
            return {true, {x0 + TARGETX, y0 + TARGETY, vx, vy, t0}};
        }
    }
    return {false, {0,0,0,0,0}};
}

// checks if a state can be created from this state which is a goal state
static pair<bool, BallState> isGoalConfigFrom(
    AngleState angles, 
    BallState ballStart, 
    AngleState& anglesp
) {
    double theta1p, theta2p, omega1p, omega2p, dx, dy, t, a1, a2;
    EffectorState effectorp;
    BallState ballt;
    pair<bool, BallState> isGoal;
    for (t = 0; t <= T_MAX; t += DT) {
        ballt = ball(ballStart, t+angles.t);
        // bias towards lower magnitude acc by starting at 0
        for (a1 = 0; a1 <= MAXACC1; a1 = a1 > 0 ? -a1 : -a1+DA) {
            // bias towards lower magnitude acc by starting at 0
            for (a2 = 0; a2 <= MAXACC2; a2 = a2 > 0 ? -a2 : -a2+DA) {
                theta1p = angles.theta1 + angles.omega1*t + 0.5*a1*t*t;
                theta2p = angles.theta2 + angles.omega2*t + 0.5*a2*t*t;
                omega1p = angles.omega1 + a1*t;
                omega2p = angles.omega2 + a2*t;
                anglesp = {theta1p, theta2p, omega1p, omega2p, t + angles.t};

                isGoal = isGoalConfig(anglesp, ballStart);
                if (isGoal.first) {
                    return isGoal;
                }
            }
        }
    }
    return {false, {0,0,0,0,0}};
}

// Euclidean distance between angle and omega configs
static double dist(AngleState start, OmegaState end) {
	double omega1 = start.omega1 - end.omega1;
    double omega2 = start.omega2 - end.omega2;
    double t = (start.t - end.t)/4; // scale time down so it doesn't overpower
	return sqrt(omega1*omega1 + omega2*omega2 + t*t);
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
static pair<int, Node*> extend(Tree* t, BallState ballStart) {
    // sample from ball start to SAMPLETIME * (time to max height) later
    OmegaState qrand = sample(ballStart.t, 
        ballStart.t + ballStart.vy / G * SAMPLETIME);
    Node *q, *qnew, *qmin = nullptr;
    double d, scale, minDist = DOUBLEMAX;
    AngleState qcheck;
    OmegaState qext;
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
    scale = EPSILON / minDist;
    if (scale >= 1) {
        qnew = new Node(arm(qmin->angles, qrand, qrand.t), qmin);
        t->V.push_back(qnew);
        return {REACHED, qnew};
    }

    qext = {
        qmin->angles.omega1 + (qrand.omega1 - qmin->angles.omega1) * scale,
        qmin->angles.omega2 + (qrand.omega2 - qmin->angles.omega2) * scale,
        qmin->angles.t + (qrand.t - qmin->angles.t) * scale
    };
    qnew = new Node(arm(qmin->angles, qext, qext.t), qmin);
    t->V.push_back(qnew);
    return {ADVANCED, qnew};
}

// RRT algorithm randomly sampling omegas and time to extend
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
		ext = extend(t, ballStart);
		s = ext.first; qnew = ext.second;
		if (s != TRAPPED) {
            // checks if can create a goal config after this state
            isGoal = isGoalConfigFrom(qnew->angles, ballStart, anglesp);
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
        {
            lock_guard<mutex> lock(worldLock);
            snapshot = world;
        }
        if (snapshot.planned > snapshot.executed * PLANAHEAD) {
            this_thread::sleep_for(chrono::milliseconds(1000));
        }
        plan = runRRT(currAngles, ballStart);
        {
            lock_guard<mutex> lock(queueLock);
            planQueue.push({plan.first, plan.second});
        }
        {
            lock_guard<mutex> lock(worldLock);
            world.planned += 1;
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
            world.executed = world.executedNext;
            if (updateBall && angleStack.empty()) {
                world.executedNext += 1;
                world.ballNext = out.second;
                updateBall = false;
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
    Text executedText, plannedText;
    Event event;
    BallState b;
    EffectorState effector;
    FloatRect bounds1, bounds2;
    
    if (!font.loadFromFile("Arial.ttf")) {
        cout << "Could not load font\n";
        return;
    }

    RenderWindow window(VideoMode(800, 600), "Arm Visualizer");
    window.setFramerateLimit(60);

    executedText.setFont(font);
    executedText.setCharacterSize(18);
    executedText.setFillColor(Color::White);

    plannedText.setFont(font);
    plannedText.setCharacterSize(18);
    plannedText.setFillColor(Color::White);

    while (window.isOpen()) {
        while (window.pollEvent(event)) {
            if (event.type == Event::Closed) {
                window.close();
            }
        }

        tnow = chrono::duration_cast<chrono::milliseconds>(
            chrono::system_clock::now().time_since_epoch()).count() 
            / 1000.0 - startTime;
        {
            lock_guard<mutex> lock(worldLock);
            world.time = tnow;
            snapshot = world;
        }
        b = ball(snapshot.ballStart, snapshot.time);

        effector = anglesToEffector(
            arm(snapshot.angleStart, snapshot.angleEnd, snapshot.time));

        Vector2f base(VISX(BASEX), VISY(BASEY));
        Vector2f joint(VISX(effector.xmid), VISY(effector.ymid));
        Vector2f end(VISX(effector.x), VISY(effector.y));

        window.clear(Color::Black);

        RectangleShape baseRect(Vector2f(50.0, 4.0f));
        baseRect.setFillColor(Color::Blue);
        baseRect.setOrigin(25.0f, 2.0f);
        baseRect.setPosition(base);

        auto link1 = makeLink(base, joint, 4.0f);
        auto link2 = makeLink(joint, end, 4.0f);
        
        window.draw(link1);
        window.draw(link2);
        window.draw(baseRect);

        effectorSize = EFFECTORRADIUS * WINDOWSCALE;
        CircleShape efffector(effectorSize);
        efffector.setFillColor(Color::Green);
        efffector.setOrigin(effectorSize, effectorSize);
        efffector.setPosition(end.x, end.y);
        window.draw(efffector);

        CircleShape circle(BALLRADIUS*WINDOWSCALE);
        circle.setFillColor(Color::Red);
        circle.setOrigin(BALLRADIUS*WINDOWSCALE, BALLRADIUS*WINDOWSCALE);
        circle.setPosition(VISX(b.x), VISY(b.y));
        window.draw(circle);

        executedText.setString("Executed Catches: " 
            + to_string(snapshot.executed));
        plannedText.setString("Planned Catches: " 
            + to_string(snapshot.planned));

        bounds1 = executedText.getLocalBounds();
        executedText.setPosition(
            window.getSize().x - bounds1.width - padding,
            padding
        );

        bounds2 = plannedText.getLocalBounds();
        plannedText.setPosition(
            window.getSize().x - bounds2.width - padding,
            padding + bounds1.height + 5
        );

        window.draw(executedText);
        window.draw(plannedText);

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
    cout << "Intial throw in 2 seconds\n";
    visualizerThread(); // run in main thread

    execution.join();
    planner.join();

    return 0;
}