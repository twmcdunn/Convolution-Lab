import java.util.ArrayList;
/**
 * Write a description of class Envelope here.
 *
 * @author (your name)
 * @version (a version number or a date)
 */
public class Envelope implements Signal
{
    double[] times, values;
    double duration;
    public Envelope(double[] t, double[] v)
    {
        times = t;
        values = v;
        duration = 1;
    }

    public void setSecondsDur(double dur){
        duration = dur;
    }

    public double getValue(double time){
        time /= duration;
        double x1 = 0;
        double x2 = 0;
        double y1 = 0;
        double y2 = 0;
        for(int i = 0; i < times.length && times[i] < time; i++){
            x1 = times[i];
            y1 = values[i];
            if(i < times.length - 1)
            {
                x2 = times[i+1];
                y2 = values[i +1];
            }
            else{
                x2 = x1 + 1;
                y2 = y1;
            }
        }
        if(x1 == x2)
            return values[0];

        double m = (y1 - y2) / (x1 - x2);
        double b = y1 - m * x1;
        return m * time + b;

    }
}
