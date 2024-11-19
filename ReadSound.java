import java.io.File;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.AudioFormat;
/**
 * Write a description of class ReadSound here.
 *
 * @author (your name)
 * @version (a version number or a date)
 */
public class ReadSound
{

    public static double[] readSoundDoubles(String fileName){
        float[] fSound = readSound(fileName);
        double[] dSound = new double[fSound.length];
        for(int i = 0; i < fSound.length; i++)
            dSound[i] = fSound[i];
        return dSound;
    }

    public static float[] readSound(String fileName){
        File  soundFile = new File(fileName);
        AudioInputStream audioInputStream = null;
        byte[] soundData = new byte[0];
        float[] normalizedData = new float[0];
        AudioFormat audioFormat = null;
        int numBytes = 0;
        try
        {
            audioInputStream = AudioSystem.getAudioInputStream(soundFile);
            audioFormat = audioInputStream.getFormat();
            int bytesPerFrame = audioInputStream.getFormat().getFrameSize();
            float frameRate = audioInputStream.getFormat().getFrameRate();
            float sampleRate = audioInputStream.getFormat().getSampleRate();
            float sampleSize = audioInputStream.getFormat().getSampleSizeInBits();
            System.out.println("isBigEndian: " + audioInputStream.getFormat().isBigEndian());
            System.out.println(audioInputStream.getFormat().getEncoding());
            System.out.println(audioInputStream.getFormat());
            if(bytesPerFrame != 2){
                System.out.println("FILE IS NOT 16-bit format");
                return null;
            }
            soundData = new byte[(int)(audioInputStream.getFrameLength() * bytesPerFrame)];
            normalizedData = new float[(int)audioInputStream.getFrameLength()];
            numBytes = audioInputStream.read(soundData);
        }
        catch (Exception ex)
        {
            System.out.println("*** Cannot process " + fileName + " ***");
            System.out.println(ex);
            System.exit(1);
        }
        int t = 0;
        float max = 0;
        int[] data = new int[normalizedData.length];

        //int[] rawIntData = new int[normalizedData.length];
        for(int i = 0; i < soundData.length; i+=2){
            /*System.out.println(soundData[i]);
            System.out.println(soundData[i + 1]);
            System.out.println();
             */

            int b1 = soundData[i];// + 128;
            if(b1 < 0){ 
                // must convert neg twos comp in little byte
                // because only big byte determines sign in 16-bit neg twos comp
                b1 = 256 + b1;
            }

            int b2 = soundData[i + 1];// + 128;
            int b2transformed = b2 << 8;
            int dataPoint = b1 | b2transformed;
            max = (float)Math.max(Math.abs(max), Math.abs(dataPoint));
            data[t] = dataPoint;
            normalizedData[t++] = (float)(dataPoint); // (double)Short.MAX_VALUE);//Short.MAX_VALUE);// - 1;

        }
        //byte[] sd2 = new byte[numBytes];
        for(int i = 0; i < normalizedData.length; i++){
            /*int c2 = data[i] / 256;
            int c1 = data[i] - c2;

            sd2[i*2] = (byte)c1;
            sd2[i*2 + 1] = (byte)(c2);
             */
            normalizedData[i] /= max;
            if(normalizedData[i] > 1)
                System.out.println("ERROR: " +normalizedData[i]);
        }

        /* SourceDataLine line = null;
        try
        {
        DataLine.Info info = new DataLine.Info(SourceDataLine.class, audioFormat);
        line = (SourceDataLine)AudioSystem.getLine(info);
        line.open(audioFormat);
        line.start();
        line.write(sd2, 0, numBytes);
        }
        catch (Exception ex)
        {
        System.out.println("*** Audio line unavailable ***");
        System.exit(1);
        }
         */
        /*
        for(int i = 0; i < sdata.length; i++){
        System.out.println(sdata[i]);
        }
         */

        return normalizedData;
    }
}
