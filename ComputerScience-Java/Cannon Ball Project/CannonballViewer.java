import javax.swing.JFrame;
import javax.swing.JOptionPane;
import java.util.Random;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import javax.swing.JComponent;

/**
   This program displays the trajectory of a cannonball.
*/
public class CannonballViewer
{
	public static final int FRAME_WIDTH = 600;
	public static final int FRAME_HEIGHT = 600;
	public static final double DELTA_T = 0.01;  //the amount of time given to the Cannonball between its displays so it can calculate
			                               //where its new position should be.

	public static final int SLEEP_DELAY = 10; //how long the program should sleep (in milliseconds) between displaying the movement of the ball
	public static Thread thread = Thread.currentThread();

	public static void main(String[] args)
	{
//		Set the terminate value to false at the begining until hit.
		boolean terminate = false;
		while(!terminate){
			JFrame frame = new JFrame();
			frame.setSize(FRAME_WIDTH, FRAME_HEIGHT);
			frame.setTitle("CannonballViewer");
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			boolean message = false;	// Shows the screen for the first time and then disappears.
			double ang;
			int target_width = 50;
			Random generator = new Random();
			int target_start = generator.nextInt(FRAME_WIDTH-target_width);

			boolean quitting = false;
		
	//		Run a while loop that finishes when it hits or i can quit it.		
			while(!quitting) {
							
				// Ask the user for the angle and velocity and angle and turn that to double type
				// for the first try it asks "Enter an angle"
				if (!message){
					String angle = JOptionPane.showInputDialog("Enter an angle:");
					ang = Double.parseDouble(angle);
					message = true;
				}else{
				// When the user misses at least once then this line is executed.
				    String miss = JOptionPane.showInputDialog("You missed. \n Try again. or type q to quit.\n Enter the initial angle");
				    if(miss.equalsIgnoreCase("Q")) {
					System.exit(0);  
				    }
				    ang = Double.parseDouble(miss);
				}
			
				String velocity = JOptionPane.showInputDialog("Enter a velocity:");
				double ivel = Double.parseDouble(velocity);

				// Creates a canball for the Cannonball class and gives the initial velocity and angle.
				 Cannonball canball = new Cannonball(ivel, ang);
				 ViewableComponents component1 = new ViewableComponents(canball, target_start, target_width);
				     
				// For the motion of the ball until it hits the floor or the target.
				 boolean motion = true;
				 while(motion) {
				     
				// Now to move the cannonball: moving on every DELTA_T set.	 
				     canball.move(DELTA_T); 
				  
				     try {  //this little snippet of code makes the program sleep for SLEEP_DELAY milliseconds
					thread.sleep(SLEEP_DELAY);
					} 
				     catch (Exception e) {
				     }
				// After awaking we add the components on the scree/ make it visible and repaint it.	
				     frame.add(component1); //Adds component1 (VisibleComponent)
				     frame.setVisible(true); //Sets frame to visible
				     component1.repaint(); // Repaints component1 to create animation

				// Checks when it landed or not.     

				     if(canball.getY()<0) 
				     motion = false;
				     else
				     motion = true;
				}
			
				// Now we check whether the CannonBall hit the target or not.
				if((canball.getX()>=target_start) && (canball.getX()<=target_start+target_width)) {
				     String congrats = JOptionPane.showInputDialog("Congratulations, you've hit the target. Play Again Y/N?");
				     if(congrats.equalsIgnoreCase("N")) {
					terminate=true; 
					quitting= true;
					
				    }
				    
				    else{
					terminate=false;
					motion=true;
					quitting =true;
					frame.remove(component1);
					frame.setVisible(false);
				    }
				}
				// If missed it will go back to the start. I have designed it so that it says the player missed and asks for the new angle.
				else {
				  
				    
				    motion=true;
				}	

			}
		}
		//Exits the program
		System.exit(0);  

	} 
}
