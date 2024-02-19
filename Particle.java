/******************************************************************************
 *  Compilation:  javac Particle.java
 *  Execution:    none
 *  Dependencies: StdDraw.java
 *      
 *  A particle moving in the unit box with a given position, velocity,
 *  radius, and mass.
 *
 ******************************************************************************/

import java.awt.Color;

/**
 *  The {@code Particle} class represents a particle moving in the unit box,
 *  with a given position, velocity, radius, and mass. Methods are provided
 *  for moving the particle and for predicting and resolvling elastic
 *  collisions with vertical walls, horizontal walls, and other particles.
 *  This data type is mutable because the position and velocity change.
 *  <p>
 *  For additional documentation, 
 *  see <a href="https://algs4.cs.princeton.edu/61event">Section 6.1</a> of 
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne. 
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class Particle {
    private static final double INFINITY = Double.POSITIVE_INFINITY;

    private double rx, ry;        // position
    private double vx, vy;        // velocity
    private int count;            // number of collisions so far
    private boolean is_brownian;  // is this the brownian particle?
    private boolean is_target;    // is this the target particle?
    private boolean is_barrier;   // is this a barrier particle?
    private boolean escaped;      // did brownian particle escaped?
    private final double radius;  // radius
    private final double mass;    // mass
    private final Color color;    // color


    /**
     * Initializes a particle with the specified position, velocity, radius, mass, and color.
     *
     * @param  rx <em>x</em>-coordinate of position
     * @param  ry <em>y</em>-coordinate of position
     * @param  vx <em>x</em>-coordinate of velocity
     * @param  vy <em>y</em>-coordinate of velocity
     * @param  radius the radius
     * @param  mass the mass
     * @param  color the color
     */
    public Particle(double rx, double ry, double vx, double vy, double radius, double mass, Color color, boolean is_brownian, boolean is_target, boolean is_barrier) {
        this.vx = vx;
        this.vy = vy;
        this.rx = rx;
        this.ry = ry;
        this.radius = radius;
        this.mass   = mass;
        this.color  = color;
        this.is_target = is_target;
        this.is_brownian = is_brownian;
        this.is_barrier = is_barrier;
        this.escaped = false; // can only be changed in the match collision case
    }
         
    /**
     * Initializes a particle with a random position and velocity.
     * The position is uniform in the unit box; the velocity in
     * either direciton is chosen uniformly at random.
     */
    public Particle() {
        rx     = StdRandom.uniform(0.0, 1.0);
        ry     = StdRandom.uniform(0.0, 1.0);
        vx     = StdRandom.uniform(-0.005, 0.005);
        vy     = StdRandom.uniform(-0.005, 0.005);
        radius = 0.02;
        mass   = 0.5;
        color  = Color.BLACK;
        is_target = false;
        is_brownian = false;
        is_barrier = false;
        escaped = false;
    }

    /**
     * Moves this particle in a straight line (based on its velocity)
     * for the specified amount of time.
     *
     * @param  dt the amount of time
     */
    public void move(double dt) {
        if(!this.is_barrier){ // so atualiza se nao for barreira
            rx += vx * dt;
            ry += vy * dt;
        }
    }

    /**
     * Draws this particle to standard draw.
     */
    public void draw() {
        StdDraw.setPenColor(color);
        StdDraw.filledCircle(rx, ry, radius);
    }

    /**
     * Returns the number of collisions involving this particle with
     * vertical walls, horizontal walls, or other particles.
     * This is equal to the number of calls to {@link #bounceOff},
     * {@link #bounceOffVerticalWall}, and
     * {@link #bounceOffHorizontalWall}.
     *
     * @return the number of collisions involving this particle with
     *         vertical walls, horizontal walls, or other particles
     */
    public int count() {
        return count;
    }

    /**
     * Returns if the escape collision involving the target exit particle
     * happened.
     */
    public boolean escaped() {
        return escaped;
    }

    /**
     * Returns if the escape collision involving the target exit particle
     * happened.
     */
    public boolean is_barrier() {
        return is_barrier;
    }

    /**
     * Returns if the particle is a target particle (exit)
     */
    public boolean is_target() {
        return is_target;
    }

    /* UNUSED
    public boolean is_micro() {
        if(!this.is_barrier && !this.is_brownian && !this.is_target) return true;
        return false;
    }
    */

    /**
     * Returns the amount of time for this particle to collide with the specified
     * particle, assuming no interening collisions.
     *
     * @param  that the other particle
     * @return the amount of time for this particle to collide with the specified
     *         particle, assuming no interening collisions; 
     *         {@code Double.POSITIVE_INFINITY} if the particles will not collide
     */
    public double timeToHit(Particle that) {
        if (this == that) return INFINITY;
        double dx  = that.rx - this.rx;
        double dy  = that.ry - this.ry;
        double dvx = that.vx - this.vx;
        double dvy = that.vy - this.vy;
        double dvdr = dx*dvx + dy*dvy;
        if (dvdr > 0) return INFINITY;
        double dvdv = dvx*dvx + dvy*dvy;
        if (dvdv == 0) return INFINITY;
        double drdr = dx*dx + dy*dy;
        double sigma = this.radius + that.radius;
        double d = (dvdr*dvdr) - dvdv * (drdr - sigma*sigma);
        //if (drdr < sigma*sigma) StdOut.println("overlapping particles");
        if (d < 0) return INFINITY;
        double prediction = -(dvdr + Math.sqrt(d)) / dvdv;
        if (prediction < 0) return INFINITY;
        return prediction;    // Foram realizadas correcoes para evitar que sejam calculados tempos negativos
    }

    /**
     * Returns the amount of time for this particle to collide with a vertical
     * wall, assuming no interening collisions.
     *
     * @return the amount of time for this particle to collide with a vertical wall,
     *         assuming no interening collisions; 
     *         {@code Double.POSITIVE_INFINITY} if the particle will not collide
     *         with a vertical wall
     */
    public double timeToHitVerticalWall() {
        double prediction;
        if      (vx > 0) prediction =  (1.0 - rx - radius) / vx;
        else if (vx < 0) prediction =  (radius - rx) / vx;
        else             return INFINITY;

        if (prediction < 0) return INFINITY;
        return prediction;    // Foram realizadas correcoes para evitar que sejam calculados tempos negativos
    }

    /**
     * Returns the amount of time for this particle to collide with a horizontal
     * wall, assuming no interening collisions.
     *
     * @return the amount of time for this particle to collide with a horizontal wall,
     *         assuming no interening collisions; 
     *         {@code Double.POSITIVE_INFINITY} if the particle will not collide
     *         with a horizontal wall
     */
    public double timeToHitHorizontalWall() {
        double prediction;
        if      (vy > 0) prediction =  (1.0 - ry - radius) / vy;
        else if (vy < 0) prediction =  (radius - ry) / vy;
        else             return INFINITY;

        if (prediction < 0) return INFINITY;
        return prediction;     // Foram realizadas correcoes para evitar que sejam calculados tempos negativos
    }

    /**
     * Updates the velocities of this particle and the specified particle according
     * to the laws of elastic collision. Assumes that the particles are colliding
     * at this instant.
     *
     * @param  that the other particle
     */
    public void bounceOff(Particle that) {
        double dx  = that.rx - this.rx;
        double dy  = that.ry - this.ry;
        double dvx = that.vx - this.vx;
        double dvy = that.vy - this.vy;
        double dvdr = dx*dvx + dy*dvy;             // dv dot dr
        double dist = this.radius + that.radius;   // distance between particle centers at collison

        // magnitude of normal force
        double magnitude = 2 * this.mass * that.mass * dvdr / ((this.mass + that.mass) * dist);

        // normal force, and in x and y directions
        double fx = magnitude * dx / dist;
        double fy = magnitude * dy / dist;

        // update velocities according to normal force
        if(!this.is_barrier){
        this.vx += fx / this.mass;
        this.vy += fy / this.mass;
        }
        if(!that.is_barrier){
        that.vx -= fx / that.mass;
        that.vy -= fy / that.mass;
        }

        // update collision counts
        if(!this.is_barrier) this.count++;
        if(!that.is_barrier) that.count++;
        if(this.is_brownian && that.is_target) this.escaped = true;
        else if(this.is_target && that.is_brownian) that.escaped = true;
    }

    /**
     * Updates the velocity of this particle upon collision with a vertical
     * wall (by reflecting the velocity in the <em>x</em>-direction).
     * Assumes that the particle is colliding with a vertical wall at this instant.
     */
    public void bounceOffVerticalWall() {
        vx = -vx;
        count++;
    }

    /**
     * Updates the velocity of this particle upon collision with a horizontal
     * wall (by reflecting the velocity in the <em>y</em>-direction).
     * Assumes that the particle is colliding with a horizontal wall at this instant.
     */
    public void bounceOffHorizontalWall() {
        vy = -vy;
        count++;
    }

    /**
     * Returns the kinetic energy of this particle.
     * The kinetic energy is given by the formula 1/2 <em>m</em> <em>v</em><sup>2</sup>,
     * where <em>m</em> is the mass of this particle and <em>v</em> is its velocity.
     *
     * @return the kinetic energy of this particle
     */
    public double kineticEnergy() {
        return 0.5 * mass * (vx*vx + vy*vy);
    }
}
