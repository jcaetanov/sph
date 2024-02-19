/******************************************************************************
 *  Compilation:  javac CollisionSystem.java
 *  Execution:    java CollisionSystem n               (n random particles)
 *                java CollisionSystem < input.txt     (from a file) 
 *  Dependencies: StdDraw.java Particle.java MinPQ.java
 *  Data files:   https://algs4.cs.princeton.edu/61event/diffusion.txt
 *                https://algs4.cs.princeton.edu/61event/diffusion2.txt
 *                https://algs4.cs.princeton.edu/61event/diffusion3.txt
 *                https://algs4.cs.princeton.edu/61event/brownian.txt
 *                https://algs4.cs.princeton.edu/61event/brownian2.txt
 *                https://algs4.cs.princeton.edu/61event/billiards5.txt
 *                https://algs4.cs.princeton.edu/61event/pendulum.txt
 *  
 *  Creates n random particles and simulates their motion according
 *  to the laws of elastic collisions.
 *
 ******************************************************************************/
import java.util.Arrays;
import java.awt.Color;

/**
 *  The {@code CollisionSystem} class represents a collection of particles
 *  moving in the unit box, according to the laws of elastic collision.
 *  This event-based simulation relies on a priority queue.
 *  <p>
 *  For additional documentation, 
 *  see <a href="https://algs4.cs.princeton.edu/61event">Section 6.1</a> of 
 *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne. 
 *
 *  @author Robert Sedgewick
 *  @author Kevin Wayne
 */
public class CollisionSystem {
    private static boolean graphic = false;

    private static final double HZ = 0.5;    // number of redraw events per clock tick

    private MinPQ<Event> pq;          // the priority queue
    private double t  = 0.0;          // simulation clock time
    private Particle[] particles;     // the array of particles

    /**
     * Initializes a system with the specified collection of particles.
     * The individual particles will be mutated during the simulation.
     *
     * @param  particles the array of particles
     */
    public CollisionSystem(Particle[] particles) {
        this.particles = particles.clone();   // defensive copy
    }

    // updates priority queue with all new events for particle a
    private void predict(Particle a, double limit) {
        if (a == null || a.is_barrier()) return; // Nao atualiza os eventos nem calcula colisoes para particulas que sao barreiras

        // particle-particle collisions
        for (int i = 0; i < particles.length; i++) {
            double dt = a.timeToHit(particles[i]);
            if (t + dt <= limit)
                pq.insert(new Event(t + dt, a, particles[i]));
        }

        // particle-wall collisions
        double dtX = a.timeToHitVerticalWall();
        double dtY = a.timeToHitHorizontalWall();
        if (t + dtX <= limit) pq.insert(new Event(t + dtX, a, null));
        if (t + dtY <= limit) pq.insert(new Event(t + dtY, null, a));
    }

    // redraw all particles
    private void redraw(double limit) {
        if(graphic) { // modo nao-grafico desenha as particulas
            StdDraw.clear();
            StdDraw.setPenColor(new Color(54, 54, 54));
            StdDraw.filledRectangle(0.5, 0.5, 300, 300);
            for (int i = 0; i < particles.length; i++) {
                particles[i].draw();
            }
            StdDraw.show();
            StdDraw.pause(20);
        }
        if (t < limit) { // mas atualiza os eventos tambem, de forma que a simulacao nao eh comprometida
            pq.insert(new Event(t + 1.0 / HZ, null, null));
        }
    }

      
    /**
     * Simulates the system of particles for the specified amount of time.
     *
     * @param  limit the amount of time
     */
    public double simulate(double limit) {
        
        // initialize PQ with collision events and redraw event
        pq = new MinPQ<Event>();
        for (int i = 0; i < particles.length; i++) {
            predict(particles[i], limit);
        }
        pq.insert(new Event(0, null, null));        // redraw event


        // the main event-driven simulation loop
        Stopwatch sw = new Stopwatch();
        while (!pq.isEmpty()) { 

            // get impending event, discard if invalidated
            Event e = pq.delMin();
            if (!e.isValid()) continue;
            Particle a = e.a;
            Particle b = e.b;

            // physical collision, so update positions, and then simulation clock
            for (int i = 0; i < particles.length; i++)
                if(!particles[i].is_barrier()) particles[i].move(e.time - t); // particulas barreira nao se movem
            t = e.time;
            StdOut.println(t);

            // process event
            if      (a != null && b != null) a.bounceOff(b);              // particle-particle collision
            else if (a != null && b == null) a.bounceOffVerticalWall();   // particle-wall collision
            else if (a == null && b != null) b.bounceOffHorizontalWall(); // particle-wall collision
            else if (a == null && b == null) redraw(limit);               // redraw event

            // update the priority queue with new collisions involving a or b
            // ONLY if the brownian particle is still inside
            if (a != null && a.escaped()) return sw.elapsedTime(); // se a particula broniana tem a flag de colisao com a saida a simulacao eh encerrada
            if (b != null && b.escaped()) return sw.elapsedTime(); // se a particula broniana tem a flag de colisao com a saida a simulacao eh encerrada
            predict(a, limit);
            predict(b, limit);
        }

        return sw.elapsedTime();
    }


   /***************************************************************************
    *  An event during a particle collision simulation. Each event contains
    *  the time at which it will occur (assuming no supervening actions)
    *  and the particles a and b involved.
    *
    *    -  a and b both null:      redraw event
    *    -  a null, b not null:     collision with vertical wall
    *    -  a not null, b null:     collision with horizontal wall
    *    -  a and b both not null:  binary collision between a and b
    *
    ***************************************************************************/
    private static class Event implements Comparable<Event> {
        private final double time;         // time that event is scheduled to occur
        private final Particle a, b;       // particles involved in event, possibly null
        private final int countA, countB;  // collision counts at event creation
                
        
        // create a new event to occur at time t involving a and b
        public Event(double t, Particle a, Particle b) {
            this.time = t;
            this.a    = a;
            this.b    = b;
            if (a != null) countA = a.count();
            else           countA = -1;
            if (b != null) countB = b.count();
            else           countB = -1;
        }

        // compare times when two events will occur
        public int compareTo(Event that) {
            return Double.compare(this.time, that.time);
        }
        
        // has any collision occurred between when event was created and now?
        public boolean isValid() {
            if (a != null && a.count() != countA) return false;
            if (b != null && b.count() != countB) return false;
            return true;
        }
   
    }


    /**
     * Unit tests the {@code CollisionSystem} data type.
     * Reads in the particle collision system from a standard input
     * (or generates {@code N} random particles if a command-line integer
     * is specified); simulates the system.
     *
     * @param args the command-line arguments
     */
    public static void main(String[] args) {

        // graphic mode
        if(args.length > 8) graphic = true;


        // number of simulations and time storage
        int run = Integer.parseInt(args[0]);
        double[] time = new double[run];
        double[] time2 = new double[run];

        // simulation proprieties
        int n1 = Integer.parseInt(args[1]);         // number of non-brownian particles
        double r1 = Double.parseDouble(args[2]);    // non-brownian radius
        double m1 = Double.parseDouble(args[3]);    // non-brownian mass (2.3258671E-26 nitrogen) (1.0E-4)
        int n2 = Integer.parseInt(args[4]);         // number of brownian particles
        double r2 = Double.parseDouble(args[5]);    // brownian and barrier initial radii
        double m2 = Double.parseDouble(args[6]);    // brownian and barrier radii (0.0025)
        double scale = Double.parseDouble(args[7]); // scale between brownian and barrier radii

        if(graphic) {
        // inicializes simulation window
        StdDraw.setCanvasSize(600, 600);

        // enable double buffering
        StdDraw.enableDoubleBuffering();
        }

        // the array of particles
        Particle[] particles;

        // create run collision systems and simulate run times
        Stopwatch total_sw = new Stopwatch();
        StdOut.print("Running simulation: ");
        for (int i = 0; i < run; i++) {
            StdOut.print((i+1)+" ");
            particles = Generator.generate(n1, r1, m1, n2, r2, m2, scale); // DEBUG needs to recieve parameters
            CollisionSystem system = new CollisionSystem(particles);
            time[i] = system.simulate(10000); // DEBUG needs to increase running time
        }

        StdOut.println("\nMean Narrow Escape Time: " + StdStats.mean(time) + " seconds.");
        StdOut.println("Total Time: " + total_sw.elapsedTime() + " seconds.");
        //StdOut.println(Arrays.toString(time)); // DEBUG
    }
      
}
