import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioSystem;
import java.io.File;
import java.io.InputStream;
import java.util.Arrays;
/**
 * Write a description of class WaveWriter here.
 *
 * @author (your name)
 * @version (a version number or a date)
 */
public class WaveWriter
{
    // instance variables - replace the example below with your own
    public StereoPcmInputStream spis;
    public float[][] df;
    public static final int SAMPLE_RATE = 48000;
    public int length = SAMPLE_RATE * 60 * 30;
    public String fileName;

    /**
     * Constructor for objects of class WaveWriter
     */
    public WaveWriter(String name)
    {
        spis = new StereoPcmInputStream();
        df = new float[2][length];
        fileName = name + ".wav";

        //render();
    }

    //render wave file
    public void render(){
        float peak = 0;
        for(int i = 0; i < df[0].length; i++){
            peak = Math.max(peak, Math.abs(df[0][i]) + Math.abs(df[1][i]));
        }
        for(int i = 0; i < df[0].length; i++){
            df[0][i] /= peak;
            df[1][i] /= peak;
        }

        for(int i = df[0].length - 1; i >= 0; i--){
            if(df[0][i] != 0 || df[1][i] != 0){
                df[0] = Arrays.copyOfRange(df[0], 0, i);
                df[1] = Arrays.copyOfRange(df[1], 0, i);
                break;
            }
        }

        spis.setDataFrames(df);

        AudioFormat af = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED, SAMPLE_RATE, 16, 2, 4,SAMPLE_RATE,false);
        try{
            AudioSystem.write(new AudioInputStream(spis, af, spis.framesToRead), AudioFileFormat.Type.WAVE, new File(fileName));
        }
        catch(Exception e){}
    }

    //basic instrument method
    public void playTone(int index, double frequency, double dynamic){
        double sec = 3;
        for(int i = index; i < length && i < index + SAMPLE_RATE * sec; i++){

            double decayFactor = ((double)(index + SAMPLE_RATE * sec - i)) / (SAMPLE_RATE * sec);
            double attackFrames = 441.0; // hundreth of second
            if(i - index < attackFrames){
                double attackFactor = (i - index) / attackFrames;
                decayFactor *= attackFactor;
            }
            df[0][i] = (float)(dynamic * decayFactor * Math.sin(frequency * Math.PI * 2 * i / SAMPLE_RATE)) + df[0][i];
        }
    }
}
