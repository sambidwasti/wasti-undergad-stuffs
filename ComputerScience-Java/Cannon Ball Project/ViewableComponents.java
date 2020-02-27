import javax.swing.JComponent;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Line2D;
import java.util.Random;
import java.awt.Color;

/**
   Draws the trajectory of a cannonball and the target.
*/
public class ViewableComponents extends JComponent
{

    private int target_start;
    private int target_width;
    
   /**
      Constructs a component that paints the flight of a cannonball and the target the ball is trying to hit
      @param ivel the initial velocity of the ball
      @param ang the angle at which the cannonball was launched
      @param t_start the starting x location of the target
      @param t_width the width of the target
   */
   public ViewableComponents(Cannonball ball, int t_start, int t_width)
   {
      this.ball = ball;
      target_start=t_start;
      target_width=t_width;
   }

   public ViewableComponents(int t_start, int t_width)
   { 
      target_start=t_start;
      target_width=t_width;
   }

   /**
    * Calls the draw() method of the CannonBall object given in the constructor.  Also draws the target on the
    * screen.
    * @param g the Graphics object on which the ball and target should be drawn.
    */
   public void paintComponent(Graphics g)
   {
      Graphics2D g2 = (Graphics2D) g;
      ball.cannon(g2);
      ball.draw(g2,this.getHeight());
      g2.draw(new Line2D.Double(target_start,this.getHeight()-1,target_start+target_width,this.getHeight()-1));
   }

   private Cannonball ball;
}
