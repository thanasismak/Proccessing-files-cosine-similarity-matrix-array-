package project2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.Arrays;
import java.lang.String;
import java.util.Scanner;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
//import org.apache.lucene.analysis.PorterStemmer;
/**
 *
 * @author Thanasis
 */

public class project2
{
    private List termsDocsArray = new ArrayList<>();//pinakas me tis lekseis apo ola ta arxeia mou meta tn epeksergasia
    private static List allTerms = new ArrayList<>(); //oles oi diaforetikes lekseis apo ola ta arxeia
    private List<double[]> tfidfDocsVector = new ArrayList<double[]>(); //mia lista me ta varoi twn leksewn pou antistoixoun ston termsdocsarray
    private List<String> names = new ArrayList<String>();
    public List author1 = new ArrayList<>(); //pinakas me ta periexomena apo to pedio author
    private final int N = 6;
    private final int M = 3;
    public String stopwards[] ={"delayed", "final", "jos", "a", "about", "above", "above", "across", "after", "afterwards", "again", "against", "all", "almost",
            "alone", "along", "already", "also", "although", "always", "am", "among", "amongst", "amoungst", "amount", "an", "and",
            "another", "any", "anyhow", "anyone", "anything", "anyway", "anywhere", "are", "around", "as", "at", "back", "be", "became",
            "because", "become", "becomes", "becoming", "been", "before", "beforehand", "behind", "being", "below", "beside", "besides",
            "between", "beyond", "bill", "both", "bottom", "but", "by", "call", "can", "cannot", "cant", "co", "con", "could", "couldnt",
            "cry", "de", "describe", "detail", "do", "done", "down", "due", "during", "each", "eg", "eight", "either", "eleven", "else",
            "elsewhere", "empty", "enough", "etc", "even", "ever", "every", "everyone", "everything", "everywhere", "except", "few",
            "fifteen", "fify", "fill", "find", "fire", "first", "five", "for", "former", "formerly", "forty", "found", "four", "from",
            "front", "full", "further", "get", "give", "go", "had", "has", "hasnt",
            "have", "he", "hence", "her", "here", "hereafter", "hereby", "herein", "hereupon", "hers", "herself",
            "him", "himself", "his", "how", "however", "hundred", "ie", "if", "in", "inc", "indeed", "interest", "into",
            "is", "it", "its", "itself", "keep", "last", "latter", "latterly", "least", "less", "ltd", "made", "many",
            "may", "me", "meanwhile", "might", "mill", "mine", "more", "moreover", "most", "mostly", "move", "much", "must",
            "my", "myself", "name", "namely", "neither", "never", "nevertheless", "next", "nine", "no", "nobody", "none",
            "noone", "nor", "not", "nothing", "now", "nowhere", "of", "off", "often", "on", "once", "one", "only", "onto",
            "or", "other", "others", "otherwise", "our", "ours", "ourselves", "out", "over", "own", "part", "per", "perhaps",
            "please", "put", "rather", "re", "same", "see", "seem", "seemed", "seeming", "seems", "serious", "several", "she",
            "should", "show", "side", "since", "sincere", "six", "sixty", "so", "some", "somehow", "someone", "something",
            "sometime", "sometimes", "somewhere", "still", "such", "system", "take", "ten", "than", "that", "the", "their",
            "them", "themselves", "then", "thence", "there", "thereafter", "thereby", "therefore", "therein", "thereupon",
            "these", "they", "thickv", "thin", "third", "this", "those", "though", "three", "through", "throughout", "thru",
            "thus", "to", "together", "too", "top", "toward", "towards", "twelve", "twenty", "two", "un", "under", "until",
            "up", "upon", "us", "very", "via", "was", "we", "well", "were", "what", "whatever", "when", "whence", "whenever",
            "where", "whereafter", "whereas", "whereby", "wherein", "whereupon", "wherever", "whether", "which", "while",
            "whither", "who", "whoever", "whole", "whom", "whose", "why", "will", "with", "within", "without", "would", "yet",
            "you", "your", "yours", "yourself", "yourselves", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "1.", "2.", "3.", "4.", "5.", "6.", "11",
            "7.", "8.", "9.", "12", "13", "14", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z",
            "terms", "CONDITIONS", "conditions", "values", "interested.", "care", "sure", ".", "!", "@", "#", "$", "%", "^", "&", "*", "(", ")", "{", "}", "[", "]", ":", ";", ",", "<", ".", ">", "/", "?", "_", "-", "+", "=",
            "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
            "contact", "grounds", "buyers", "tried", "said,", "plan", "value", "principle.", "forces", "sent:", "is,", "was", "like",
            "discussion", "tmus", "diffrent.", "layout", "area.", "thanks", "thankyou", "hello", "bye", "rise", "fell", "fall", "psqft.", "http://", "km", "miles","\\\""};

    public void parsefile(String file) throws IOException, FileNotFoundException{


       File[] allfiles = new File(file).listFiles();
    BufferedReader input = null;
       
       for(File f: allfiles){ 
           if(f.getName().endsWith(".bib")){
              input=new BufferedReader(new FileReader(f));
              StringBuilder sb = new StringBuilder();
    String s = null;
    names.add(f.getName().replaceAll(".bib", ""));
              while((s = input.readLine())!=null){
                 sb.append(s);
              }
String[] tokenizedTitles = sb.toString().split(" title     =");
String[] tokenizedJournal = sb.toString().split(" journal   =");
String[] tokenizedBooktitles = sb.toString().split(" booktitle =");
String[] tokenizedAuthor = sb.toString().split("  author    = ");//author
String tempArray = "";
String tempArray1 = "";
             for(int i = 1; i<tokenizedTitles.length; i++){
                tempArray=tempArray+(tokenizedTitles[i].substring(tokenizedTitles[i].indexOf("{")+1, tokenizedTitles[i].indexOf("},")))+" ";//kratame tis lekseis eksrontas tis parenthesis
             }
             for(int i = 1; i<tokenizedBooktitles.length; i++){
                tempArray=tempArray+(tokenizedBooktitles[i].substring(tokenizedBooktitles[i].indexOf("{")+1, tokenizedBooktitles[i].indexOf("},")))+" ";
             }
             for(int i = 1; i<tokenizedJournal.length; i++){
                tempArray=tempArray+(tokenizedJournal[i].substring(tokenizedJournal[i].indexOf("{")+1, tokenizedJournal[i].indexOf("},")))+" ";
             }
             for(int i = 1; i<tokenizedAuthor.length; i++){//author
                 tempArray1=tempArray1+(tokenizedAuthor[i].substring(tokenizedAuthor[i].indexOf("{")+1, tokenizedAuthor[i].indexOf("},")))+" ";
             }
            String[] allterms2 = tempArray.toLowerCase().replaceAll(",|[{(.)}^]|:|/", " ").replaceAll("-", " ").replaceAll("\\s+", " ").split(" ");
String[] allterms4 = tempArray1.toLowerCase().replaceAll(",|[{(.)}^]|:|/", " ").replaceAll("-", " ").replaceAll("\\s+", " ").split(" "); //author
ArrayList<String> allterms5 = new ArrayList<>();
            for(String term: allterms4){//author
                boolean flag = true;
               for(String term2:stopwards){
                   if(term2.equals(term)){
                       flag=false;
                       break;
                    }
                }
                if(flag){
                    allterms5.add(term);
                }
            }
            author1.add(allterms5); //author
            PrintWriter output5 = new PrintWriter("author.txt", "UTF-8");
//System.out.println(author1);
output5.println(author1);
            output5.close();
            ArrayList<String> allterms3 = new ArrayList<>();
            for(String term: allterms2){
                boolean flag = true;
               for(String term2:stopwards){
                   if(term2.equals(term)){
                       flag=false;
                       break;
                    }
                }
                if(flag){
                    allterms3.add(term);
                }
            }
            termsDocsArray.add(allterms3);
            for(String term: allterms3){
                if(!allTerms.contains(term)){
                    allTerms.add(term);
                }
            }
           }
        }
       //System.out.println(allTerms.size());//megethos sunoliko
    }
    public double idfCalculator(List allTerms, String termToCheck)
{
    double count = 0;
    for (int i = 0; i < allTerms.size(); i++)
    {
        ArrayList<String> temparray = new ArrayList<>();
        temparray = (ArrayList<String>)allTerms.get(i);
        for (int j = 0; j < temparray.size(); j++)
        {

            if (temparray.get(j).equals(termToCheck))
            {
                count++;
                break; //mas endiafferei na uparxei i leksi sto sugkekrimeno arxeio!
            }
        }
    }
    return Math.log(allTerms.size() / count);
}
public double tfCalculator(ArrayList<String> totalterms, String termToCheck)
{
    double count = 0.0;  //to count the overall occurrence of the term termToCheck
    for (String s : totalterms)
    {
        if (s.equals(termToCheck))
        {
            count++;
        }
    }
    return Math.log(1 + count / totalterms.size());
}
public void TfIdfCalculator()
{
    double tf = 0, idf = 0, tfidf = 0;
    boolean flag = false;
    for (int i = 0; i < termsDocsArray.size(); i++)
    {
        double[] tfidfvectors = new double[allTerms.size()];
        ArrayList<String> temparray = new ArrayList<>();
        temparray = (ArrayList<String>)termsDocsArray.get(i);
        for (int j = 0; j < allTerms.size(); j++)
        {
            tf = tfCalculator(temparray, allTerms.get(j).toString());
            idf = idfCalculator(termsDocsArray, allTerms.get(j).toString());
            tfidf = tf * idf;
            tfidfvectors[j] = tfidf;
        }
        tfidfDocsVector.add(tfidfvectors);
    }
}
//ERWTIMA1 : FIND MOST IMPORTANT WORDS FOR EVERY TEACHER
public void erwtima1() throws FileNotFoundException, UnsupportedEncodingException{
        PrintWriter output = new PrintWriter("prof-description.txt", "UTF-8");
        for(int i = 0; i<names.size(); i++){
            double[] tempterms = new double[allTerms.size()];
double[] finalterms = new double[N];
int[] ind = new int[N];
            for (int ii = 0; ii<N; ii++){ind[ii] = N+1;}
            String s1 = "";
s1=s1+names.get(i)+" ";
            tempterms=Arrays.copyOfRange(tfidfDocsVector.get(i), 0, tfidfDocsVector.get(i).length);
            Arrays.sort(tempterms);
            finalterms=Arrays.copyOfRange(tempterms, tempterms.length-N, tempterms.length);
            int index = N;
            for (int l = 0; l<=index; l++){
                double temp = finalterms[index - 1];
finalterms[index - 1] = finalterms[l];
                finalterms[l] = temp;
                index--;
            }
            for(int j = 0; j<N; j++){
                for(int k = 0; k<tempterms.length; k++){
                    boolean flag = false;
                    if(finalterms[j]==tfidfDocsVector.get(i)[k]){
                        for(int l = 0; l<N; l++){
                            if (ind[l]==k){ flag= true;}
                        }
                        if (!flag){
                        ind[j]=k;
                        break;   
                        }
                    }
                }
                s1=s1+allTerms.get(ind[j])+","+String.format("%.3f",finalterms[j])+ " ";
            }
            output.println(s1);
            
        }
        output.close();
    }
    public double CosineSimilarity(double[] docVector1, double[] docVector2)
{
    double dotProduct = 0.0;
    double magnitude1 = 0.0;
    double magnitude2 = 0.0;
    double cosineSimilarity = 0.0;

    for (int i = 0; i < docVector1.length; i++)
    {
        dotProduct += docVector1[i] * docVector2[i];  //a.b
        magnitude1 += Math.pow(docVector1[i], 2);  //(a^2)
        magnitude2 += Math.pow(docVector2[i], 2); //(b^2)
    }

    magnitude1 = Math.sqrt(magnitude1);//sqrt(a^2)
    magnitude2 = Math.sqrt(magnitude2);//sqrt(b^2)

    if (magnitude1 != 0.0 | magnitude2 != 0.0)
    {
        cosineSimilarity = dotProduct / (magnitude1 * magnitude2);
    }
    else
    {
        return 0.0;
    }
    return cosineSimilarity;
}
//ERWTISI2: FIND TEACHER BY A KEY QUESTION
public void erwtima2() throws IOException
{
    BufferedReader s1=new BufferedReader(new InputStreamReader(System.in));
        String query = s1.readLine();
String[] querytmp = query.toLowerCase().replaceAll(",|[{(.)}^]|:|/", " ").replaceAll("-", " ").replaceAll("\\s+", " ").split(" ");
ArrayList<String> querytmp2 = new ArrayList<>();
        for(String term: querytmp){
               boolean flag = true;
               for(String term2:stopwards){
                   if(term2.equals(term)){
                       flag=false;
                       break;
                    }
                }
                if(flag){
                        if(allTerms.contains(term)){
                            querytmp2.add(term);
                        }
                }
        }
        double tf = 0, idf = 0, tfidf = 0;
double[] tfidfvectors = new double[allTerms.size()];
        for(int j = 0; j<allTerms.size(); j++){
            tf=tfCalculator(querytmp2, allTerms.get(j).toString());
            idf=1; //log(N/n) exoume 1 erwtisi=N exoume toulaxiston 1 emfanisi=n se alli periptwsi xrisimopoiontas ton tupou tis diafaneias to varos mas tha evgene 0 ara k 0 to tfidf 
            tfidf=tf* idf;
tfidfvectors[j]=tfidf;
        }
        double[] querysim = new double[names.size()];
int count = 0;
        for(double[] vector:tfidfDocsVector){
            querysim[count]=(CosineSimilarity(tfidfvectors, vector));
            count++;
        }
        double[] tempterms = new double[names.size()];
tempterms=Arrays.copyOfRange(querysim, 0,querysim.length);
        Arrays.sort(tempterms);
        int[] ind2 = new int[names.size()];
        for (int i = 0; i<names.size(); i++){ind2[i] = names.size()+1;}
        int index = names.size();
        for (int l = 0; l<=index; l++){
            double temp = tempterms[index - 1];
tempterms[index - 1] = tempterms[l];
            tempterms[l] = temp;
            index--;
        }
        for(int j = 0; j<tempterms.length; j++){
            boolean flag = false;
            for(int i = 0; i<querysim.length; i++){
                if(querysim[i]==tempterms[j]){
                    for(int l = 0; l<names.size(); l++){
                        if (ind2[l]==i){ flag= true;}
                    }
                    if(!flag){
                        ind2[j]=i;
                    }
                }
            }
        }
        String fname = "results";
        for(String s:querytmp2){
           fname=fname+"-"+s;
        }
        PrintWriter output2 = new PrintWriter(fname + ".txt", "UTF-8");
        for(int k = 0; k<names.size(); k++){
           String s2 = "";
s2=s2+names.get(ind2[k])+","+String.format("%.3f",tempterms[k]);
            output2.println(s2);
        }
        output2.close();
    }
    //ERWTIMA3 : FIND SIMILAR TEACHER BY KEY WORD
    public void erwtima3() throws UnsupportedEncodingException, FileNotFoundException{
        PrintWriter output3 = new PrintWriter("similar_profs.txt", "UTF-8");
List<double[]> simMatrix = new ArrayList<double[]>(); 
        for(int i = 0; i<names.size(); i++){
            double[] arraytmp = new double[names.size()];
            for(int j = 0; j<names.size(); j++){
                arraytmp[j]=CosineSimilarity(tfidfDocsVector.get(i),tfidfDocsVector.get(j));//ftiaxnw to simMatrxi (tin kathe grammi)
            }
            simMatrix.add(arraytmp);
        }
        for(int i = 0; i<names.size(); i++){
            double[] arraytmp2 = new double[names.size()];
arraytmp2=Arrays.copyOf(simMatrix.get(i), names.size());//tin kathe grami tou simMatrix
            Arrays.sort(arraytmp2);
            double[] arraytmp3 = new double[M];
arraytmp3=Arrays.copyOfRange(arraytmp2,arraytmp2.length-M-1,arraytmp2.length-1);//sorted array me M times pio omoies
            int index = M;
            for (int l = 0; l<=index; l++){ //antistrefw tis theseis tou pinaka mou apo megalutero-mikrotero!
                double temp = arraytmp3[index - 1];
arraytmp3[index - 1] = arraytmp3[l];
                arraytmp3[l] = temp;
                index--;
            }
            int[] ind = new int[M];
            for (int ii = 0; ii<M; ii++){ind[ii] = M+1;}
            String s1 = "";
s1=s1+names.get(i)+" ";
            for(int j = 0; j<M; j++){
                for(int k = 0; k<names.size(); k++){
                    boolean flag = false;
                    if(arraytmp3[j]==simMatrix.get(i)[k]){
                        for(int l = 0; l<M; l++){
                            if (ind[l]==k){ flag= true;}
                    }
                        if (!flag){
                            ind[j]=k;
                        break;   
                        }
                    }
                }
                s1=s1+names.get(ind[j])+","+String.format("%.3f", arraytmp3[j] )+ " ";
            }
            output3.println(s1);
        }
        output3.close();
    }
    //na valw ena orio sts epanalupseis
    public void kmeans()throws FileNotFoundException, UnsupportedEncodingException{ //8-12 epanalupseis sunithws!
       final int L = allTerms.size() / 2;
PrintWriter outputY = new PrintWriter("tety.txt", "UTF-8");
PrintWriter outputT = new PrintWriter("tett.txt", "UTF-8");
double[] seedY = new double[L]; //kentro A
double[] seedT = new double[L]; //kentro B (krataei suntetagmenes kapoiou simeiou giauto einai double[]
Set<Integer> clusterT = new HashSet<Integer>(); //lista apo integer
Set<Integer> clusterY = new HashSet<Integer>(); //lista apo integer
List<double[]> tfidfL = new ArrayList<double[]>(); //lista 26xL me ta varh twn leksewn
        //kratame ta L stoixeia apo tfidfvectors
        for(int i = 0; i<26; i++){
            double[] temptfidfL = tfidfDocsVector.get(i);
double[] temptfidfL2 = new double[L];
            for(int j = 0; j<L; j++){
                temptfidfL2[j]=temptfidfL[j];
            }
            tfidfL.add(temptfidfL2);
        }
        //pernoume 1 kathigiti apo tilepoikinwnies kai 1 apo upologistes
        seedY=Arrays.copyOfRange(tfidfL.get(24), 0, tfidfL.get(24).length); //pernoume ton prwto kathigiti ws kentro gia upologistes "athanasiadou"
        seedT=Arrays.copyOfRange(tfidfL.get(2), 0, tfidfL.get(2).length);// pernoume ton deutero kathigiti ws kentro gia tis tilepoikinwnies "blionas"
        clusterY.add(24); //WALLACE
        clusterT.add(2); //boucouvalas
        //k-means
        boolean stopflag = false;
Set<Integer> clusterTprev = new HashSet<Integer>();
Set<Integer> clusterYprev = new HashSet<Integer>();
        while(!stopflag){
            clusterTprev.clear();
            clusterYprev.clear();
            clusterTprev.addAll(clusterT); //krataw tous kathigites kathe cluster/omadas
            clusterYprev.addAll(clusterY);
            //vima1
            clusterT.clear();
            clusterY.clear();
            for(int i = 0; i<26; i++){
                double distY = eudist(tfidfL.get(i), seedY); // apostasi anamesa stis suntetagmenes kathe kathigiti kai tou kentrou
double distT = eudist(tfidfL.get(i), seedT);
                if(distY>= distT){
                    clusterT.add(i);
                }else
                    clusterY.add(i);
            }
            //vima2
            double[] sumY = new double[L];
double[] sumT = new double[L];
                for(int j = 0; j<sumY.length; j++){
                    sumY[j]=0;
                    sumT[j]=0;
                }
                for(int i:clusterT){
                    sumT=sumarrays(sumT, tfidfL.get(i));
                }
                for(int i = 0; i<seedT.length; i++){
                    seedT[i]= sumT[i] / clusterT.size();
                }
                for(int i: clusterY){
                    sumY=sumarrays(sumY, tfidfL.get(i));
                }
                for(int j = 0; j<seedY.length; j++){
                    seedY[j]=sumY[j] / clusterY.size();
                }
                String sy = "";
String st = "";
                for(int i: clusterT){
                    st=st+names.get(i)+","+String.format("%.3f",eudist(tfidfL.get(i),seedT))+ " ";
                }
                for(int j: clusterY){
                    sy=sy+names.get(j)+","+String.format("%.3f",eudist(tfidfL.get(j),seedY))+ " "; 
                }
                if((clusterT.equals(clusterTprev))&& (clusterY.equals(clusterYprev))){ //an einai idia ta cluster(omades) me tin proigoumeni epanalipsu -an einai idia ta clusterT tha einai kai ta clusterY 
                    stopflag=true;
                }
                 outputT.println(st);
                 outputY.println(sy);
        }
        outputY.close();
        outputT.close();
    }
    public double[] sumarrays(double x[], double y[])
{
    double[] c = new double[x.length];
    for (int i = 0; i < x.length; i++)
    {
        c[i] = x[i] + y[i];
    }
    return c;
}
public double eudist(double x[], double y[])
{
    double Sum = 0;
    double distance = 0;
    for (int i = 0; i < x.length; i++)
    {
        Sum = Sum + pow((x[i] - y[i]), 2.0);
        distance = sqrt(Sum);
    }
    return distance;
}
public void askisi2erw2()
{
    ArrayList<double[]> B = new ArrayList<double[]>();
    for (int i = 0; i < 26; i++)
    { //arxikopoihsh pinaka B
        for (int j = 0; j < 26; j++)
        {
            B.get(i)[j] = 0;
        }
    }
    for (int i = 0; i < 26; i++)
    { //kathigites
        int count = 0;
        ArrayList<String> temp = new ArrayList<>();
        temp.add((String)author1.get(i)); //pernoume to 1o string px athanasiadou
        for (int j = 0; j < temp.size(); j++)
        {//lekseis(georgia)
         //System.out.println("eeesa2");
            String current = temp.get(j); //kathe leksei apo to temp
            for (int ii = 0; ii < 26; ii++)
            {// kathigites
             //System.out.println("eeesa");
                for (int k = 0; k < names.size(); k++)
                { //elegxei tin kathe leksi me ta onomata pou exoume
                  //System.out.println("eeesa1");
                    if (current.equals(names.get(k)))
                    { //
                        if (k != i)
                        { //diaforo apo ton kathigiti pou elegxw twra
                            B.get(i)[k] += 1; //+1
                            count += 1; //to athroisma twn varwn pou uparxoun stin grammi
                        }
                    }
                }

            }
        }
        for (int jj = 0; jj < 26; jj++)
        {
            B.get(i)[jj] = B.get(i)[jj] / count;
        }
    }
    System.out.println(B);
}


public static void main(String args[])throws FileNotFoundException, IOException, NullPointerException
    {
      project2 skt = new project2();
skt.parsefile("papers");
      skt.TfIdfCalculator();
      skt.erwtima1();
      //skt.erwtima2();
      //skt.erwtima3();
      skt.kmeans();
      skt.askisi2erw2();
    }
}    

