import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import java.awt.Polygon;
import java.awt.Color;
/**
   This class simulates a cannonball fired at an angle.
*/
public class Cannonball
{
   /**
      Constructs a Cannonball 
      @param anIvel the initial velocity of the ball
      @param anAng the angle at which the cannonball was launched
   */
   public Cannonball(double ivel, double ang)
   {
      angle = ang;
      posX = 0;
      posY = 0;
      velX = Math.cos(Math.toRadians(ang)) * ivel;
      velY = Math.sin(Math.toRadians(ang)) * ivel;
   }
   
   /**
    * Draws the Cannonball on the graphics object.  The height of the window is
    * also passed in so that the ball's correct y location can be drawn.
    * @param graphics the Graphics2D object on which this object draws itself
    * @param height the heigth of the viewable window, used to calculate the ball's
    * y location.
    */
   public void draw(Graphics2D graphics,int height) {
	 final double RADIUS = 5;
         double x = getX(); 
         double y = height - getY();
	if (getY()<=0)	graphics.setColor(Color.MAGENTA);	// For some extra special effect.
	else	graphics.setColor(Color.RED);	
         Ellipse2D.Double circle = new Ellipse2D.Double(
            x - RADIUS, y - RADIUS, 2 * RADIUS, 2 * RADIUS);
         graphics.fill(circle);
   }

   /**
      Updates the position and velocity of this cannon ball 
      after a given time interval.
      @param deltaT the time interval
   */
   public void move(double deltaT)
   {
      posX = posX + velX * deltaT;
      posY = posY + velY * deltaT;
      velY = velY - G * deltaT;
   }


   /**
	Draw a Cannon and point it according to the angle
   */
   public void cannon(Graphics2D g2)
   {
	int size = 15;	
	int x=(int)(size*Math.cos((Math.PI*getAngle()/180)));
        int y=(CannonballViewer.FRAME_HEIGHT-30)-(int)(size*Math.sin((Math.PI*getAngle()/180)));
	Polygon cannonpoint = new Polygon();
	cannonpoint.addPoint(5,CannonballViewer.FRAME_HEIGHT-30);
	cannonpoint.addPoint(0,CannonballViewer.FRAME_HEIGHT-35);
	cannonpoint.addPoint(x,y-5);
	cannonpoint.addPoint(x+5,y);
	if (getY()<=0)	g2.setColor(Color.BLACK);		// for some extra color flashy effect.
	else	g2.setColor(Color.GREEN);	
	g2.fill(cannonpoint);
	Ellipse2D.Double circle = new Ellipse2D.Double(0-10, CannonballViewer.FRAME_HEIGHT-30-10, 20, 20);
	g2.fill(circle);
   }
  /**
      Gets the x position of this cannon ball.
      @return the horizontal position
   */
   public double getX()
   {
      return posX;
   }

   /**
      Gets the y position of this cannon ball.
      @return the vertical position
   */
   public double getY()
   {
      return posY;
   }
   /** 
      Gets the angle
      @ return the angle
   */
   public double getAngle()
   {
       return angle;
   }

   private static final double G = 9.81;
   
   private double posX;
   private double posY;
   private double velX;
   private double velY;
   private double angle;

}
