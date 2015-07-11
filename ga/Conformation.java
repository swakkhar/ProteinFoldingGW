/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package ga;

/**
 *
 * @author Administrator
 */
import java.util.ArrayList;
public class Conformation  implements Comparable{
    boolean mutate;
    boolean pullsinglepoint;
    boolean cornerflip;
    boolean pull;
    boolean tilt;
    boolean crankshaft;

            private int size;
            public AminoAcid[] aminoAcid;
	    public int[] fitness={0,0};
            double energy=0.0;
	    public String encodedString;
            //public Operation[] operation;
            Conformation()
	    {
	        size=Global.seqLength;
                this.aminoAcid=null;
	        //this.fitness=0;
                this.encodedString="";
                //this.operation=new Operation[size];
	    }
            Conformation(AminoAcid[] AA,int Fitness[])
	    {
                size=Global.seqLength;
                aminoAcid=new AminoAcid[size];
                int i=0;
                //int[] op={-1,-1,-1,-1,-1};
	        for (AminoAcid a:AA)
	        {
	            this.aminoAcid[i]=new AminoAcid(a.x,a.y,a.z,a.operation);
	            i++;

	        }
	        this.fitness[0]=Fitness[0];
                this.fitness[1]=Fitness[1];
                this.energy=(double)this.fitness[0]/1000.0;
                this.encodedString="";



	    }
	    Conformation(AminoAcid[] AA,int[] Fitness,String EncodedSequence)
	    {
                size=Global.seqLength;
                aminoAcid=new AminoAcid[size];
                int i=0;
	        for (AminoAcid a:AA)
	        {
	            this.aminoAcid[i]=new AminoAcid(a.x,a.y,a.z,a.operation);
	            i++;
	        }
	        this.fitness[0]=Fitness[0];
                this.fitness[1]=Fitness[1];
                this.energy=(double)this.fitness[0]/1000.0;
                this.encodedString=EncodedSequence;
	    }

            Conformation(Conformation _c)
	    {
                size=Global.seqLength;
                aminoAcid=new AminoAcid[size];
                int i=0;

	        for (AminoAcid a:_c.aminoAcid)
	        {
	            this.aminoAcid[i]=new AminoAcid(a.x,a.y,a.z,a.operation);
	            i++;
	        }
	        //this.fitness=_c.fitness;
                this.fitness[0]=_c.fitness[0];
                this.fitness[1]=_c.fitness[1];
                this.energy=(double)this.fitness[0]/1000.0;
                this.encodedString=_c.encodedString;
	    }


	    public int getMaxX()
	    {
	    int X = this.aminoAcid[0].x;
	        for (int i = 1; i < size; i++)
	        {
	            if (this.aminoAcid[i].x >X)
	            {
	                X = this.aminoAcid[i].x;
	            }
	        }
	        return X;
	    }
	    public int getMinX()
	    {
	        int X = this.aminoAcid[0].x;
	        for (int i = 1; i < size; i++)
	        {
	            if (this.aminoAcid[i].x <X)
	            {
	                X = this.aminoAcid[i].x;
	            }
	        }
	        return X;
	    }
	    public int getMaxY()
	    {
	        int Y = this.aminoAcid[0].y;
	        for (int i = 1; i < size; i++)
	        {
	            if (this.aminoAcid[i].y >Y)
	            {
	                Y = this.aminoAcid[i].y;
	            }
	        }
	        return Y;
	    }
	    public int getMinY()
	    {
	        int Y = this.aminoAcid[0].y;
	        for (int i = 1; i < size; i++)
	        {
	            if (this.aminoAcid[i].y <Y)
	            {
	                Y = this.aminoAcid[i].y;
	            }
	        }
	        return Y;
	    }
	    public int getMaxZ()
	    {
	        int Z = this.aminoAcid[0].z;
	        for (int i = 1; i < size; i++)
	        {
	            if (this.aminoAcid[i].z >Z)
	            {
	                Z = this.aminoAcid[i].z;
	            }
	        }
	        return Z;
	    }
	    public int getMinZ()
	    {
	        int Z = this.aminoAcid[0].z;
	        for (int i = 1; i < size; i++)
	        {
	            if (this.aminoAcid[i].x <Z)
	            {
	                Z = this.aminoAcid[i].z;
	            }
	        }
	        return Z;
	    }
    public Conformation randConformation(){
        boolean valid = false;
        String sequence = "";
        int[] fitness={0,0};
        AminoAcid[] aa = new AminoAcid[size];
        while (!valid) {
            int x = 0;
            int y = 0;
            int z = 0;

            aa[0] = new AminoAcid(x, y, z);
            ArrayList<Node3d> occupiedNode = new ArrayList<Node3d>();
            occupiedNode.add(new Node3d(x, y, z));

            for (int i = 1; i < size; i++) {
                int p = 0;
                ArrayList<Node3d> al = new ArrayList<Node3d>();
                al = Global.vecArrays.get(p);
                int k = 0;
                for (int TRY = 0; TRY < 12; TRY++) {
                    k =Common.getRandomNext(0,12);
                    Node3d node = al.get(k);
                    x = aa[i - 1].x + node.x;
                    y = aa[i - 1].y + node.y;
                    z = aa[i - 1].z + node.z;
                    valid = true;
                    for (Node3d el : occupiedNode) {
                        if (el.x == x && el.y == y && el.z == z) {
                            valid = false;
                            break;
                        }
                    }
                    if (valid) {
                        break;
                    }
                }
                if (!valid) {
                    break;
                }
                aa[i] = new AminoAcid(x, y, z);
                occupiedNode.add(new Node3d(x, y, z));
            }
        }

        RelativeZone rz=new RelativeZone();
        sequence=rz.encode(aa);
        //int fitnesss=Common.calculateFitness(aa);
        Common.calculateFitness(aa,fitness);
        Conformation c = new Conformation(aa, fitness, sequence);
        return c;
    }
public Conformation initialConformation()
{
    int confSize=Global.seqLength-1;
        int _i=0,_ii;
	int ra;
	int r_size=Common.getRandomNext(0,4)+1;//% (/*confSize/15*/4)+1;
        String str="";
        char[] conf=new char[confSize];
	//cout<<r_size<<endl;
	while(_i<confSize)
	{
		ra= (Common.getRandomNext(1,r_size));
		for(_ii=0;_ii<ra && _i<confSize;_ii++)
		{
//			conf[_i++]='a'+0;
                        conf[_i++]='a'+0;

		}
		if(_i>=confSize) break;
		//conf[_i++]='a'+1;
                conf[_i++]='a'+1;
		if(_i>=confSize) break;
		//ra= (rand()%r_size)+1;
                ra= (Common.getRandomNext(1,r_size));
		for(_ii=0;_ii<ra && _i<confSize;_ii++)
		{
			//conf[_i++]='a'+2;
                        conf[_i++]='a'+9;
		}
		if(_i>=confSize) break;



		//conf[_i++]='a'+4;
                conf[_i++]='a'+1;

		if(_i>=confSize) break;


		//ra= (rand()%r_size)+1;
                ra= (Common.getRandomNext(1,r_size));
		for(_ii=0;_ii<ra && _i<confSize;_ii++)
		{
			//conf[_i++]='a'+2;
                        conf[_i++]='a'+0;
		}
		if(_i>=confSize) break;
//		conf[_i++]='a'+3;
                conf[_i++]='a'+1;
		if(_i>=confSize) break;
		//ra= (rand()%r_size)+1;
                ra= (Common.getRandomNext(1,r_size));
		for(_ii=0;_ii<ra && _i<confSize;_ii++)
		{
//			conf[_i++]='a'+0;
                        conf[_i++]='a'+9;
		}
		if(_i>=confSize) break;


		//conf[_i++]='a'+4;
                conf[_i++]='a'+1;


		if(_i>=confSize) break;
	}
        for(char cr:conf){
            str+=Character.toString(cr);
        }
//        System.out.println(str.length());
//        System.out.println(str);
        RelativeZone rz=new RelativeZone();
        Conformation c = rz.decode(str);
        return c;

}
    public int compareTo(Object obj) {
            if (this.fitness[1] == ((Conformation)obj).fitness[1])
                return 0;
            else if (this.fitness[1] > ((Conformation)obj).fitness[1])
                return 1;
            else
                return -1;

    }
}