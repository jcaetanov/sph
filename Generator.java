/******************************************************************************
 *  Compilation:  javac Generator.java
 *  Execution:    java Generator n1 r1 m1 n2 r2 m2 scale
 *                [n1] number of non-browninan particles
 *                [r1] radius of non-browninan particles
 *                [m1]   mass of non-browninan particles
 *                [n2] number of browninan particles
 *                [r2] radius of browninan particles
 *                [m2]   mass of browninan particles
 *             [scale] scale between brownian and barrier radii
 *            
 *  Dependencies: StdRandom.java 
 *  
 *	Based on user given molecular properties, generates initial conditions
 *	for a Brownian Motion simulation
 *
 ******************************************************************************/
import java.util.Arrays;
import java.awt.Color;
import java.lang.Math;

public class Generator {
	private static final double kb = 1.38064852E-23; // Boltzmann constant [m2 kg s^-2 K^-1]

	// debug: allows particles a breathing space so collision errors do not occur
	private static double breather() {
		return StdRandom.uniform(0.0000000000773721, 0.0000000002910635); // experimental: may need tweaking
	}
	
	private static double[] maxwell(int n, double temp, double mass){
		// velocities array
		double[] vel = new double[n];

		// velocity component population
		for (int i = 0; i < n; i++) 
			vel[i] = Math.sqrt((kb * temp) / mass) * StdRandom.gaussian(0, 1);
		
		//StdOut.println(Arrays.toString(vel)); // DEBUG
		return vel;
	}
	
	public static Particle[] generate(int n1, double r1, double m1, int n2, double r2, double m2, double scale) {

		// barrier proprieties (includes exit)
		double r3 = r2 * scale;
		double d3 = r3 * 2;
		int n3 = (int) Math.floor(1 / d3);

		// particles array
		Particle[] particles = new Particle[(n1+n2+n3)]; // n1 = non-brownian; n2 = brownian; n3 = barrier.

		// particle proprieties
		double rx;
        double ry;
        double vx;
        double vy;
        double radius;
        double mass;
        Color  color;

        // system proprieties
        double energy = 0.00000125; // kelvin

        // color palletes
        Color[] warm  = {new Color(229,  33,  35), new Color(231, 51,  34), new Color(233, 81,  32), new Color(235, 109,  33), new Color(235, 136,  35), new Color(238, 163,  35)}; //n1_color
        Color[] cold  = {new Color( 44,   0, 212), new Color( 61,  9, 213), new Color( 85, 22, 214), new Color(116,  24, 218), new Color(126,  23, 219), new Color(145,  24, 220),
                         new Color(172,  26, 222), new Color(202, 26, 226), new Color(227, 30, 224)}; //n3_color
        Color   exit  =  new Color(19, 189, 59);
        int pallete; // iterator
		
		// generates n1 non-brownian particles: uniform position maxwell velocity
		radius = r1;
		mass   = m1;

		pallete = 0; // iterates through colors pallete
		for(int i = 0; i < n1; i++) {
			rx = StdRandom.uniform(r1, 1 - (d3 + r3)); // generates random position from [d3, 1-d3)
			ry = StdRandom.uniform(r1, 1 - r3); // concentrates particles inside barrier, avoiding collision
			double[] maxwell_vx = maxwell(n1, energy, mass); // USING NITROGEN MASS 
			double[] maxwell_vy = maxwell(n1, energy, mass);
			particles[i] = new Particle(rx, ry, maxwell_vx[i], maxwell_vy[i], radius, mass, warm[pallete], false, false, false); // non-brownian particle
			if(pallete == warm.length-1) pallete = 0; // color iteration
			else pallete++; // color iteration
		}

		// generates n2 brownian particles
		rx     = 0.5; // center position
		ry     = 0.5; 
		vx     = 0.0; // still
		vy     = 0.0;
		radius = r2;
		mass   = m2;
		color  = new Color(33, 180, 250); //n2_color
		
		for(int j = n1; j < (n1+n2); j++)
			particles[j] = new Particle(rx, ry, vx, vy, radius, mass, color, true, false, false); // brownian particle
	
		// generates barrier particles
		vx     = 0.0; // still
		vy     = 0.0;
		radius = r3;
		mass   = 100000.0; // still

		ry = radius; // y iterator
		pallete = 0; // color iterator
		int exit_pivot = (int) Math.floor(((n1+n2+n3) - (n1+n2)) / 2); // sets exit particle in the average center

		for(int k = (n1+n2); k < (n1+n2+n3); k++) { // left
			if(k == (n1+n2)+exit_pivot)
				particles[k] = new Particle((1-radius), ry, vx, vy, (radius-breather()), mass, exit, false, true, true); // exit particle
			else {
				particles[k] = new Particle((1-radius), ry, vx, vy, (radius-breather()), mass, cold[pallete], false, false, true); // barrier particles
				if(pallete == cold.length-1) pallete = 0; // color iteration
				else pallete++; // color iteration
			}
			ry += d3 + breather();
		}

		return particles;
	}

	/* test client*/
	public static void main(String[] args) {
		// argument reading
		int n1 = Integer.parseInt(args[0]);
		double m1 = Integer.parseInt(args[1]);
		double r1 = Integer.parseInt(args[2]);
		int n2 = Integer.parseInt(args[3]);
		double m2 = Integer.parseInt(args[4]);
		double r2 = Integer.parseInt(args[5]);
		double scale = Integer.parseInt(args[6]);

		// array generator
		Particle[] particles = generate(n1, r1, m1, n2, r2, m2, scale);

		// debug
		StdOut.println(Arrays.toString(particles));
	}

}