import java.util.*; //using * pulls in all classes from package
import javax.swing.JFrame;

/**Animated viewer for my Clock.  This class contains only the main method to
 * run and view the animated clock.
 * Created by Louis Oliphant
 * Date 1/29/2010
 */

public class ClockViewer {

    /**Initial size of the square viewing window.
     * 
     */
    public static final int WINSIZE=600;
    
    /**The amount of time (in milliseconds) that the program sleeps between updates
     * of the graphics window.
     */
    public static final int DELAY=100;

    
    /**Runs the ClockViewer.  Displays an animated clock which constantly
     * updates itself with the current time.
     */
    public static void main(String args[]) {
        
        //create the frame and set its title and closing behaviour
        JFrame frame = new JFrame();
        frame.setSize(WINSIZE,WINSIZE);
        frame.setTitle("My Clock");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        //create a calendar object and pull out the current time
        //(in hours, minutes, and seconds
        Calendar calendar=new GregorianCalendar();
        int hour=calendar.get(Calendar.HOUR);
        int minute=calendar.get(Calendar.MINUTE);
        int second=calendar.get(Calendar.SECOND);
	int hourdigital=calendar.get(Calendar.HOUR_OF_DAY);
        
        //create my Clock object, place it in the frame, and make frame visible
        Clock clock = new Clock(hour,minute,second,hourdigital);
        frame.add(clock);
        frame.setVisible(true);
        
        //The thread holds access to this running program
        Thread thread=Thread.currentThread();

        //keep looping forever (until the user clicks close)
        while(true) {
            //grab the current time
            calendar = new GregorianCalendar();
            hour=calendar.get(Calendar.HOUR);
            minute=calendar.get(Calendar.MINUTE);
            second=calendar.get(Calendar.SECOND);
            hourdigital=calendar.get(Calendar.HOUR_OF_DAY);

            //tell the clock the current time
            clock.setTime(hour,minute,second,hourdigital);
            
            //tell the clock to repaint itself
            clock.repaint();
            
            //now sleep the running program for a while
            try {
                thread.sleep(DELAY);
            } catch(Exception e) {
            }
        }
    }
    
}
