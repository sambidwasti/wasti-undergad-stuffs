import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Ellipse2D;
import javax.swing.JComponent;
import java.awt.Color;
import java.awt.Rectangle;			


/**A graphical representation of an analog clock, with a clock face and hands for the hour,
 * minute, and second.  The clock resizes itself to keep itself centered in the window
 * and at full screen size.
 */
class Clock extends JComponent {

    private int xCenter; //the x coordinate of the center point of the clock
    private int yCenter; //the y coordinate of the center point of the clock
    private int radius;  //the radius of the outer circle (the clockface size)

    private int hour;   	//holds the hour setting for the hour hand
    private int minute; 	//holds the minute setting for the minute hand
    private int second; 	//holds the second setting for the second hand
    private int hourdigital; 	//holds the hour for the 24 hr clock in digital

    /**Constructor for the Clock class.  Initializes the time of the clock to the hour,
     * minute, and second.
     * @param hr the hour setting for the hour hand
     * @param m the minute setting for the minute hand
     * @param s the second setting for the second hand
     */
    public Clock(int h, int m, int s, int dh) {
        setTime(h,m,s,dh);
    }
    
    /**Updates the clocks time.  Sets the instance fields to be the same
     * as those passed in as parameters.
     * @param h the hour setting for the hour hand
     * @param m the minute setting for the minute hand
     * @param s the second setting for the second hand
     */
    public void setTime(int h, int m, int s, int dh) {
        this.hour=h;
        this.minute=m;
        this.second=s;
	this.hourdigital=dh;
    }
    
    /**Draws the clock on the given graphics object.  This class extends the JComponent class.
     * This method is called (eventually) when the repaint() method is called.
     * @param graphics the graphics object on which this Clock object should draw itself.
     */
    public void paintComponent(Graphics graphics) {
	

        
        //find the midpoint of the viewable screen (xCenter,yCenter) and
        //set clock radius to be half of the smaller dimension
        int width=this.getWidth();
        int height=this.getHeight();
        if (width<height) {
            radius=width/2;
        } else {
            radius=height/2;
        }
        xCenter=width/2;
        yCenter=height/2;
        
        Color clrclock = Color.BLACK;			// Define the color of the Clock
        Graphics2D g=(Graphics2D)graphics;
        
        //draw the circle around the middle point
        Ellipse2D.Double circle=new Ellipse2D.Double(xCenter-radius,yCenter-radius,radius*2,radius*2);
	g.setColor(clrclock);
        g.fill(circle);
        
        //draw digital readout version of clock
        drawDigitalReadout(g);

        //draw the numbers around the edge
        drawNumbers(g);
        
        //draw the hourHand
        drawHourHand(g);
        
        //draw the minuteHand
        drawMinuteHand(g);
        
        //draw the secondHand
        drawSecondHand(g);
        
    }
    
    
    private void drawDigitalReadout(Graphics2D g) {
	int xcord = xCenter-42;	
	int ycord = yCenter+50;	
	Rectangle box1 = new Rectangle(xcord-3,ycord-3,91,26);
	Rectangle box2 = new Rectangle(xcord, ycord, 85, 20);	// define a box so that the display is there.
	g.setColor(Color.MAGENTA);
	g.draw(box1);
	g.draw(box2);

	String digitaltext = hourdigital+" : "+minute+" : "+second;
	g.drawString(digitaltext, xcord+10,  ycord+15);	
    
    }
    
    private void drawHourHand(Graphics2D g) {
        //set the radius for how long to draw the hand based upon the radius of the outer circle
        int hrHandRadius=(int)(radius*0.85*2/3);
        
        //figure out the end point of the hour hand
        //recall that a full circle is 360 degrees
        //in radians this is 2*pi
        //and the min hand changes every min, i.e. 1/12 of a rotation for every hour
        int xLocation=xCenter+(int)(hrHandRadius*Math.sin((2*Math.PI*hour/12.0)+(2*Math.PI*minute/(12*60.0))));
        int yLocation=yCenter-(int)(hrHandRadius*Math.cos((2*Math.PI*hour/12.0)+(2*Math.PI*minute/(12*60.0))));
        
        //draw the hour hand
	g.setColor(Color.RED);
        g.drawLine(xCenter,yCenter,xLocation,yLocation);
    }
    
    
    private void drawMinuteHand(Graphics2D g) {
        //set the radius for how long to draw the hand based upon the radius of the outer circle
        int minHandRadius=(int)(radius*0.80);
        
        //figure out the end point of the min hand
        //recall that a full circle is 360 degrees
        //in radians this is 2*pi
        //and the min hand changes every min, i.e. 1/60 of a rotation for every min 
        int xLocation=xCenter+(int)(minHandRadius*Math.sin(2*Math.PI*minute/60.0));
        int yLocation=yCenter-(int)(minHandRadius*Math.cos(2*Math.PI*minute/60.0));
        
        //draw the minute hand
	g.setColor(Color.YELLOW);
        g.drawLine(xCenter,yCenter,xLocation,yLocation);
    }
    
    
    /**A helper function which draws the second hand on the clock face.
     * @param g the Graphics2D object on which to draw the second hand.
     */
    private void drawSecondHand(Graphics2D g) {
        
        //set the radius for how long to draw the hand based upon the radius of the outer circle
        int secHandRadius=(int)(radius*0.85);
        
        //figure out the end point of the second hand
        //recall that a full circle is 360 degrees
        //in radians this is 2*pi
        //and the second hand changes every second, i.e. 1/60 of a rotation for every second
        int xLocation=xCenter+(int)(secHandRadius*Math.sin(2*Math.PI*second/60.0));
        int yLocation=yCenter-(int)(secHandRadius*Math.cos(2*Math.PI*second/60.0));
        
        //draw the second hand
	g.setColor(Color.BLUE);
        g.drawLine(xCenter,yCenter,xLocation,yLocation);
    }
    
    
    /**A helper function that draws the numbers around the clock face.
     * @param g the Graphics2D object on which the numbers are placed.
     */
    private void drawNumbers(Graphics2D g) {
        
        //set radius for where to place the numbers based upon the radius of the outer circle
        int numRadius=(int)(radius*0.95);
        
        //height and width of a single number (these are just a guesstimate)
        int charHeight=5;
        int charWidth=3;
        
        //now loop around the clock face and place the numbers
        for(int hourNum=1;hourNum<=12;hourNum++) {
            
            //recall that a full circle has 360 degrees.
            //in radians this is 2*pi
            //and the numbers need to be placed around the circle every 1/12 of the way
            int xLocation=xCenter+(int)(numRadius*Math.sin(2.0*Math.PI*hourNum/12.0));
            int yLocation=yCenter-(int)(numRadius*Math.cos(2.0*Math.PI*hourNum/12.0));
            
	    g.setColor(Color.CYAN);
	    if (hourNum < 10){
		Ellipse2D.Double cirhourNum=new Ellipse2D.Double(xLocation-10,yLocation-10,20,20);
		g.draw(cirhourNum);
	    }	
	    else{
		Ellipse2D.Double cir1hourNum=new Ellipse2D.Double(xLocation-7,yLocation-10,22,20);
		g.draw(cir1hourNum);
	    }	
            //draw the number and the appropriate location (offset by its width and height)
	   
	    
            g.drawString(""+hourNum,xLocation-charWidth,yLocation+charHeight);
        }
    }
}
