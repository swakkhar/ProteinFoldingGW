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
import java.util.Random;
//import java.util.Enumeration;
import java.io.*;
public class Global {
    public static int energydistselected=1;
    public static int energydistributions[][][]=new int[2][280][280];

    public static String aaseq;
    public static int energyfunction[][]=new int[256][256];    

    public static Random ranDOM;//=new Random(290001);
    //public static Random random;//=new Random(290001);
    public static long rndSeed=0;
    public static boolean SomethingWrong = false;
         
    public static String pid="";
    

    public static boolean traceOn=false;
    public static boolean logging=false;
        
    public static FileOutputStream fout;
    public static FileOutputStream foutlog;
    public static FileOutputStream foutoplog;
    public static String inputFile="";
    //=============================================================================
    public static int populationSize=20;
    public static int reservedSize=200;
    public static long solExplored=0;

    public static int maxTry = 10;
    public static long startTimeInSec = 0;

    public static double bestFitness;

    public static int seqLength;
    public static int restartCount=0;
    public static int maxRestart=0;

    public static int generationCount = 0;
    public static int totalGenerationFitness = 0;
    
    public static String twelveDirections="abcdefghijkl";

    public static int xoverAttempts=0;
    public static int xoverSucceeds=0;
    public static int mutationAttempts=0;
    public static int mutationSucceeds=0;
    public static int pullAttempts=0;
    public static int pullSucceeds=0;
    public static int flipAttempts=0;
    public static int flipSucceeds=0;
    public static int tiltAttempts=0;
    public static int tiltSucceeds=0;
    
    public static int discardedTwins=0;

    public static int greeddivision[]={100,100,500,1000,1500,2000,3477};
    public static int greeddivsize=7;

    public static int greedparam[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    public static int seqgreedlevel;

    public static int iterationparam[]={160,160,150,125,118,105,85,55,39,33,21,16,12,11,7,5,4};//{100,100,92,80,75,66,55,38,28,30,18,12};
    //{130,130,120,110,102,90,84,54,37,31,22,15};//{100,100,92,80,75,66,55,38,28,30,18,12};

    public static int stagparam[]={50,50,50,50,45,45,40,35,30,25,20,15,15,15,15,15,15};

    public static int seqstagparam;
    public static int distswitchparam;

    public static String[][] benchmarkprimary={
        {"4RXN", "mkkytctvcgyiynpedgdpdngvnpgtdfkdipddwvcplcgvgkdqfeevee"},
        {"1ENH", "rprtafsseqlarlkrefnenrylterrrqqlsselglneaqikiwfqnkraki"},
        {"4PTI", "rpdfcleppytgpckariiryfynakaglcqtfvyggcrakrnnfksaedcmrtcgga"},
        {"2IGD", "mtpavttyklvingktlkgetttkavdaetaekafkqyandngvdgvwtyddatktftvte"},
        {"1YPA", "mktewpelvgkavaaakkvilqdkpeaqiivlpvgtivtmeyridrvrlfvdkldniaqvprvg"},
        {"1R69", "sissrvkskriqlglnqaelaqkvgttqqsieqlengktkrprflpelasalgvsvdwllngtsdsnvr"},
        {"1CTF", "aaeektefdvilkaagankvavikavrgatglglkeakdlvesapaalkegvskddaealkkaleeagaevevk"},
        {"3MX7", "mtdlvavwdvalsdgvhkiefehgttsgkrvvyvdgkeeirkewmfklvgketfyvgaaktkatinidaisgfayeytleingkslkkym"},
        {"3NBM", "snaskelkvlvlcagsgtsaqlanaineganltevrviansgaygahydimgvydliilapqvrsyyremkvdaerlgiqivatrgmeyihltkspskalqfvlehyq"},
        {"3MQO", "paidyktafhlapiglvlsrdrviedcndelaaifrcaradligrsfevlypssdeferigerispvmiahgsyaddrimkraggelfwchvtgraldrtaplaagvwtfedlsatrrva"},
        {"3MRO", "snalsaseerfqlavsgasaglwdwnpktgamylsphfkkimgyedhelpdeitghresihpddrarvlaalkahlehrdtydveyrvrtrsgdfrwiqsrgqalwnsagepyrmvgwimdvtdrkrdedalrvsreelrrl"},
        {"3PNX", "gmenkkmnlllfsgdydkalasliianaaremeievtifcafwgllllrdpekasqedkslyeqafssltpreaeelplskmnlggigkkmllemmkeekapklsdllsgarkkevkfyacqlsveimgfkkeelfpevqimdvkeylknalesdlqlfi"},
        {"3MSE", "ispnvlnnmksymkhsnirniiinimahelsvinnhikyinelfykldtnhngslshreiytvlasvgikkwdinrilqaldindrgnitytefmagcyrwkniestflkaafnkidkdedgyisksdivslvhdkvldnndidnfflsvhsikkgiprehiinkisfqefkdymlstf"},
        {"3MR7", "snaerrlcailaadmagysrlmernetdvlnrqklyrrelidpaiaqaggqivkttgdgmlarfdtaqaalrcaleiqqamqqreedtprkeriqyriginigdivledgdifgdavnvaarleaisepgaicvsdivhqitqdrvsepftdlglqkvknitrpirvwqwvpdadrdqshdpqpshvqh"},
        {"3NO6", "mtfskelreasrpiiddiyndgfiqdllagklsnqavrqylradasylkeftniyamlipkmssmedvkflveqiefmlegeveahevladfinepyeeivkekvwppsgdhyikhmyfnafarenaaftiaamapcpyvyavigkramedpklnkesvtskwfqfystemdelvdvfdqlmdrltkhcsetekkeikenflqstiherhffnmayinekweyggnnne"},
        {"3NO3", "mnlkstlllllclmmagmvaakdntkviahrgywktegsaqnsirsleraseigaygsefdvhltadnvlvvyhdndiqgkhiqsctydelkdlqlsngeklptleqylkrakklknirlifelkshdtpernrdaarlsvqmvkrmklakrtdyisfnmdackefirlcpksevsylngelspmelkelgftgldyhykvlqshpdwvkdckvlgmtsnvwtvddpklmeemidmgvdfittdlpeetqkilhsraq"},
        {"3ON7", "mkletidyraadsakrfveslretgfgvlsnhpidkelveriytewqaffnseaknefmfnrethdgffpasisetakghtvkdikeyyhvypwgripdslranilayyekantlasellewietyspdeikakfsiplpemianshktllrilhyppmtgdeemgairaaahedinlitvlptanepglqvkakdgswldvpsdfgniiinigdmlqeasdgyfpstshrvinpegtdktksrislplflhphpsvvlserytadsylmerlrelgvl"},
     };

    public static ArrayList<Node3d> vecA=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecB=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecC=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecD=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecE=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecF=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecG=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecH=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecI=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecJ=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecK=new ArrayList<Node3d>();
    public static ArrayList<Node3d> vecL=new ArrayList<Node3d>();

    public static ArrayList<ArrayList> vecArrays=new ArrayList<ArrayList>();

    //public static int minX = 0, minY = 0, maxX = 0, maxY = 0;
    public static int TRY = 0;
    public static int reTry = 0;
    public static int fitnessImprove = 0;
    public  static int[][][] nextMove;


    public static char[] baseDir=new char[12];
    public static int[][] baseVec;
    public static int[][][] baseMat;
    public static int [][][] invrMat;
    public static Conformation bestFitConf;
    

    public  static Conformation[] arrayInitialPopulation=new Conformation[populationSize];
    public  static Conformation[] arrayCurrentPopulation=new Conformation[populationSize];

    public static ArrayList<Conformation> alReservedSol=new ArrayList<Conformation>();


    public enum operatorType{
        // 0=crossover,1=classic mutation,2=pullmove,3=corner-flip,4=tilt
        crossOver(0),classicMutation(1),pullmove(2),cornerFlip(3),tiltMove(4);
        private final int type;
        operatorType(int type) {
            this.type=type;
        }
        public int opType(){
            return this.type;
        }
    }


    public static class NeighborData {
        char dir;
        Node3d vec;
        Node3d[] mat = new Node3d[3];
        Node3d[] invmat = new Node3d[3];
        NeighborData(char dir,Node3d vec,Node3d[] baseMat,Node3d[] invMat){
            this.dir=dir;
            this.vec=vec;
            this.mat=baseMat;
            this.invmat=invMat;
        }
    }
    
    public static Node3d[] vector3Dfcc=new Node3d[12];

    public static ArrayList<Node3d> occupiedNode=new ArrayList<Node3d>();
}
