
import java.io.*;
import java.text.DateFormat;
import java.util.*;
import java.util.concurrent.*;


import modi.*;
import moda.OneMOD;
import moda.MultiMOD;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.JDOMException;
import org.jdom.input.SAXBuilder;

import moda.DPHeap;
import msutil.IsobaricTag;
import msutil.MSMass;
import msutil.PGraph;
import msutil.ProtCutter;
import processedDB.*;
import scaniter.MSMScan;
import scaniter.ScanIterator;

public class MODPlus {
    private static final String[] message = {
            "[Error] Cannot read any MS/MS scan from input dataset.\r\n" +
                    "[Error] Check consistency between input file and its format.",

            "[Error] Cannot read any protein from input database.\r\n" +
                    "[Error] Check input fasta format.",

            "[Error] One fixed modification per amino acid can be allowed.\r\n" +
                    "[Error] Check specfied fixed modifications.",

            "[Error] Unsupported character set in your search parameter",

            "[Error] Required field is empty.\r\n" +
                    "[Error] Required fields : MS/MS Data, Database",

            "[Error] Wrong usage.\r\n" +
                    "[Error] Re-confirm it.",

            "[Error] Not defined"
    };
    static boolean dynamicPMCorrection = false, multiBlind = true;
    static int numHeatedPeptides = 50;

    public static void main(String[] args) throws Exception {
        Constants.engine = "modplus";
        Constants.engineVersion = "hyu";


        int availableCores = Runtime.getRuntime().availableProcessors();


        System.out.println("************************************************************************************");
        System.out.println("Modplus (version " + Constants.engineVersion + ") - Identification of post-translational modifications");
        System.out.println("Release Date: 2025");
        System.out.println("Available CPU Cores: " + availableCores);
        System.out.println("************************************************************************************");
        System.out.println();


        run(args[0]);
    }


    protected static int set_parameter(String Prixparam) throws Exception {

        System.out.println("Reading parameter.....");
        System.out.println("DEBUG: cpu-AvailableProcessors: " + Runtime.getRuntime().availableProcessors());
        System.out.println();

        Document doc;
        try {
            doc = new SAXBuilder().build(Prixparam);
        } catch (JDOMException e) {
            System.out.println(message[3]);
            return 1;
        } catch (IOException e) {
            System.out.println(message[5]);
            return 5;
        }

        Element search = doc.getRootElement();
        Constants.runDate = DateFormat.getDateInstance().format(new Date());
        if (search.getAttributeValue("user") != null) {
            Constants.runUser = search.getAttributeValue("user");
        }
        if (search.getAttributeValue("title") != null) {
            Constants.runTitle = search.getAttributeValue("title");
        } else Constants.runTitle = String.valueOf(System.currentTimeMillis());

        Element dataset = search.getChild("dataset");
        if (dataset != null) {
            Constants.SPECTRUM_LOCAL_PATH = dataset.getAttributeValue("local_path");
            if (Constants.SPECTRUM_LOCAL_PATH == "") {
                System.out.println(message[4]);
                return 4;
            }

            String type = dataset.getAttributeValue("format");
            if (type.compareToIgnoreCase("mgf") == 0) Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MGF;
            else if (type.compareToIgnoreCase("pkl") == 0) Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.PKL;
            else if (type.compareToIgnoreCase("ms2") == 0) Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MS2;
            else if (type.compareToIgnoreCase("dta") == 0) Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.DTA;
            else if (type.compareToIgnoreCase("mzxml") == 0)
                Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.MZXML;
            else if (type.compareToIgnoreCase("zip") == 0)
                Constants.SPECTRA_FILE_TYPE = Constants.spectra_format.ZIPDTA;

            Constants.INSTRUMENT_NAME = dataset.getAttributeValue("instrument");
            if (Constants.INSTRUMENT_NAME.equals("QTOF")) Constants.INSTRUMENT_TYPE = Constants.msms_type.QTOF;
            else Constants.INSTRUMENT_TYPE = Constants.msms_type.TRAP;
        }
        System.out.print("Input datasest : " + Constants.SPECTRUM_LOCAL_PATH);
        System.out.println(" (" + Constants.SPECTRA_FILE_TYPE + " type)");

        Element database = search.getChild("database");
        if (database != null) {
            Constants.PROTEIN_DB_LOCAL_PATH = database.getAttributeValue("local_path");
            if (Constants.PROTEIN_DB_LOCAL_PATH == "") {
                System.out.println(message[4]);
                return 4;
            }
        }
        System.out.println("Input database : " + Constants.PROTEIN_DB_LOCAL_PATH);

        Element enzyme = search.getChild("enzyme");//DEPRECATED
        if (enzyme != null) {
            String enzymeName = enzyme.getAttributeValue("name");
            String cut = enzyme.getAttributeValue("cut");
            String sence = enzyme.getAttributeValue("sence");
            Constants.protease = new ProtCutter(enzymeName, cut, sence);
        }//*/

        Element com_enzyme = search.getChild("combined_enzyme");
        if (com_enzyme != null) {
            String enzymeName = com_enzyme.getAttributeValue("name");
            String nn = com_enzyme.getAttributeValue("nterm_cleave");
            String cc = com_enzyme.getAttributeValue("cterm_cleave");
            Constants.protease = new ProtCutter(enzymeName, nn, cc, true);
        }

        Constants.alkylatedToCys = 0; //DEPRECATED
        Element cys_alkylated = search.getChild("cys_alkylated");
        if (cys_alkylated != null) {
            Constants.alkylationMethod = cys_alkylated.getAttributeValue("name");
            Constants.alkylatedToCys = Double.parseDouble(cys_alkylated.getAttributeValue("massdiff"));
            AminoAcid.modifiedAminoAcidMass('C', Constants.alkylatedToCys);
            MSMass.modifiedAminoAcidMass('C', Constants.alkylatedToCys);
        }

        Element instrument_resolution = search.getChild("instrument_resolution");
        if (instrument_resolution != null) {
            Constants.MSResolution = ("high".compareToIgnoreCase(instrument_resolution.getAttributeValue("ms")) == 0) ? 1 : 0;
            if (Constants.MSResolution == 1) System.out.println("High resolution MS!!");
            Constants.MSMSResolution = ("high".compareToIgnoreCase(instrument_resolution.getAttributeValue("msms")) == 0) ? 1 : 0;
            if (Constants.MSMSResolution == 1) System.out.println("High resolution MS/MS!!");
        }

        Element parameters = search.getChild("parameters");
        if (parameters != null) {
            Element param;
            param = parameters.getChild("enzyme_constraint");
            if (param != null) {
                Constants.missCleavages = Integer.parseInt(param.getAttributeValue("max_miss_cleavages"));
                Constants.numberOfEnzymaticTermini = Integer.parseInt(param.getAttributeValue("min_number_termini"));
                if (Constants.numberOfEnzymaticTermini > 2) Constants.numberOfEnzymaticTermini = 2;
            }

            param = parameters.getChild("isotope_error");
            if (param != null) {
                if (param.getAttributeValue("min_C13_number") != null)
                    Constants.minNoOfC13 = Integer.parseInt(param.getAttributeValue("min_C13_number"));

                if (param.getAttributeValue("max_C13_number") != null)
                    Constants.maxNoOfC13 = Integer.parseInt(param.getAttributeValue("max_C13_number"));

                if (Constants.maxNoOfC13 == 0 && param.getAttributeValue("increment_per_dalton") != null)
                    Constants.rangeForIsotopeIncrement = Integer.parseInt(param.getAttributeValue("increment_per_dalton"));
            }

            param = parameters.getChild("peptide_mass_tol");
            if (param != null) {
                if (param.getAttributeValue("unit").compareToIgnoreCase("ppm") == 0) {
                    Constants.PPMTolerance = Double.parseDouble(param.getAttributeValue("value"));
                } else {
                    Constants.precursorTolerance = Constants.precursorAccuracy = Double.parseDouble(param.getAttributeValue("value"));
                }
            }

            param = parameters.getChild("fragment_ion_tol");
            if (param != null) {
                Constants.fragmentTolerance = Double.parseDouble(param.getAttributeValue("value"));
            }
            param = parameters.getChild("modified_mass_range");
            if (param != null) {
                Constants.minModifiedMass = Double.parseDouble(param.getAttributeValue("min_value"));
                Constants.maxModifiedMass = Double.parseDouble(param.getAttributeValue("max_value"));
            }
        }

        Element protocol = search.getChild("protocol");
        if (protocol != null) {
            System.out.print("Protocol Description: ");
            Element isobaric = protocol.getChild("isobaric_labeling");
            if (isobaric != null) {
                if (isobaric.getAttributeValue("reagent") != null) {
                    Constants.isobaricTag = isobaric.getAttributeValue("reagent");
                    Constants.reporterMassOfIsobaricTag = IsobaricTag.getReporterMasses(isobaric.getAttributeValue("reagent"));
                    if ("".compareTo(Constants.isobaricTag) != 0)
                        System.out.print(Constants.isobaricTag + " Labelled" + ((Constants.reporterMassOfIsobaricTag == null) ? " (NOT Supported)" : " (Supported)"));
                }
            }
            Element modEnrich = protocol.getChild("modification_enrichment");
            if (modEnrich != null) {
                if (modEnrich.getAttributeValue("mod") != null) {
                    Constants.enrichedModification = modEnrich.getAttributeValue("mod");
                    if ("".compareTo(Constants.enrichedModification) != 0)
                        System.out.print(" & " + Constants.enrichedModification + " Enriched" +
                                (("Acetyl".compareToIgnoreCase(Constants.enrichedModification) == 0 || "Phospho".compareToIgnoreCase(Constants.enrichedModification) == 0) ? " (Supported)" : " (NOT Supported)"));
                }
            }
            System.out.println();
        }

        Element modifications = search.getChild("modifications");
        Constants.variableModifications = new PTMDB();
        Constants.fixedModifications = new PTMDB();

        if (modifications != null) {
            double[] fixedAA = new double[26];

            Element fixed = modifications.getChild("fixed");
            if (fixed != null) {
                if (Constants.fixedModifications.setFixedModificatinos(fixed, fixedAA) == 0) {
                    System.out.println(message[2]);
                    return 2;
                }
            }
            if (!Constants.fixedModifications.isEmpty())
                System.out.println("Fixed modifications : " + Constants.fixedModifications.size() + " selected");

            Element variable = modifications.getChild("variable");
            if (variable != null) {
                Constants.PTM_FILE_NAME = variable.getAttributeValue("local_path");
                boolean canBeModifiedOnFixedAA = variable.getAttributeValue("canBeModifiedOnFixedAA").equals("1");
                Constants.canBeModifiedOnFixedAA = canBeModifiedOnFixedAA;
                if (Constants.PTM_FILE_NAME != null) {
                    Constants.variableModifications.setVariableModificatinos(Constants.PTM_FILE_NAME, fixedAA, canBeModifiedOnFixedAA);
                }
                Constants.variableModifications.setVariableModificatinos(variable, fixedAA, canBeModifiedOnFixedAA);

                if (canBeModifiedOnFixedAA) {
                    for (PTM p : Constants.fixedModifications) {
                        Constants.variableModifications.add(
                                new PTM(Constants.variableModifications.size(), "De-" + p.getName(), "",
                                        -p.getMassDifference(), 0, p.getResidue(), p.getPTMPosition(), (p.getAbbAA() == 'C') ? 1 : 0));
                    }
                }
                if (variable.getAttributeValue("multi_mods") != null && variable.getAttributeValue("multi_mods").equals("0")) {
                    Constants.maxPTMPerGap = Constants.maxPTMPerPeptide = 1;
                }
            }
            if (!Constants.variableModifications.isEmpty()) {
                System.out.print("Variable modifications : " + Constants.variableModifications.size() + " selected (");
                Constants.variableModifications.setPTMDiagnosticIon();
                if (Constants.maxPTMPerPeptide == 1) System.out.println("one modification per peptide)");
                else System.out.println("multiple modifications per peptide)");
            }
        }
        Constants.variableModifications.constructPTMLookupTable();

        Element decoy_search = search.getChild("decoy_search");
        if (decoy_search != null) {
            if (decoy_search.getAttributeValue("checked") != null) {
                if ("1".compareTo(decoy_search.getAttributeValue("checked")) == 0) {
                    Constants.targetDecoy = 1;
                    System.out.println("Decoy search checked");
                }
            }
        }

        Element multistages_search = search.getChild("multistages_search");
        if (multistages_search != null) {
            if (multistages_search.getAttributeValue("checked") != null) {
                if ("1".compareTo(multistages_search.getAttributeValue("checked")) == 0) {
                    Constants.firstSearchProgram = multistages_search.getAttributeValue("program");
                    System.out.println("MultiStages Search checked " + Constants.firstSearchProgram);
                }
            }
        }

        Element mod_map = search.getChild("mod_map");
        if (mod_map != null) {
            if (mod_map.getAttributeValue("checked") != null) {
                if ("1".compareTo(mod_map.getAttributeValue("checked")) == 0) {
                    System.out.println("MODMap checked");
                }
            }
        }

        Constants.adjustParameters();

        System.out.println();
        return 0;
    }


    public static void run(String arg) throws Exception {


        try {
            if (set_parameter(arg) != 0) return;
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        try {
            File analPath = new File(Constants.SPECTRUM_LOCAL_PATH);
            if (analPath.isDirectory()) {
                String type = Constants.SPECTRA_FILE_TYPE.toString().toLowerCase();
                for (File file : analPath.listFiles()) {
                    if (file.getName().endsWith(type)) {
                        Constants.SPECTRUM_LOCAL_PATH = file.getPath();
                        System.out.println("Input dataset: " + Constants.SPECTRUM_LOCAL_PATH);
                        modplus_mod_search();
                    }
                }
                System.out.println("End of process");
            } else {
                modplus_mod_search();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    static int modplus_mod_search() throws Exception {
        System.out.println("Starting MODPlus for modification search!");

        Constants.MAX_TAG_SIZE = 100;
        Constants.minTagLength = 2;
        Constants.minTagLengthPeptideShouldContain = 3;
        Constants.tagChainPruningRate = 0.4;

        String identifier = Constants.SPECTRUM_LOCAL_PATH;
        identifier = identifier.substring(0, identifier.lastIndexOf('.'));


        ScanIterator scaniter = ScanIterator.get(Constants.SPECTRUM_LOCAL_PATH, Constants.SPECTRA_FILE_TYPE);
        if (scaniter == null || scaniter.size() == 0) {
            System.out.println("Failed to read msms spectra file");
            return 1;
        }
        System.out.println(scaniter.size() + " scans");


        System.out.print("Reading protein database.....  ");
        StemTagTrie ixPDB = new StemTagTrie(Constants.PROTEIN_DB_LOCAL_PATH);
        if (ixPDB.getSizeOfEntries() == 0) {
            System.out.println("Failed to read protein fasta file");
            return 1;
        }
        System.out.println();

        int NUM_THREADS = Runtime.getRuntime().availableProcessors();

        ThreadPoolExecutor executor = new ThreadPoolExecutor(
                NUM_THREADS,                 // 코어 스레드 수
                NUM_THREADS,                 // 최대 스레드 수
                0L,                          // 유휴 스레드 유지 시간
                TimeUnit.MILLISECONDS,
                new LinkedBlockingQueue<>(18_000), // 작업 큐 (크기 제한)
                Executors.defaultThreadFactory(),
                new ThreadPoolExecutor.CallerRunsPolicy() // 큐가 가득 찼을 때 호출 스레드에서 실행
        );

        long startTime = System.currentTimeMillis();


        boolean considerIsotopeErr = (Constants.maxNoOfC13 != 0 || Constants.precursorTolerance > 0.50001) ? true : false;

        int totalSize = scaniter.size();

        // Create a list to hold all scan
        List<ArrayList<MSMScan>> allScans = new ArrayList<>();
        while (scaniter.hasNext()) {
            allScans.add(scaniter.getNext());
        }

        // Create a structure to hold results in order
        ConcurrentHashMap<Integer, ResultEntry> results = new ConcurrentHashMap<>();


        // task submit per thread
        CountDownLatch latch = new CountDownLatch(allScans.size());
        for (int i = 0; i < allScans.size(); i++) {  // outer for loop
            final int scanIndex = i;
            executor.execute(() -> {
                try {
                    int selected = -1;
                    ArrayList<AnsPeptide> candidates = null;
                    ArrayList<AnsPeptide> tp;

                    final ArrayList<MSMScan> scanTarget = allScans.get(scanIndex);
                    final int sz = scanTarget.size();
                    for (int j = 0; j < sz; j++) {
                        final int JJ = j;
                        Spectrum spectrum = scanTarget.get(JJ).getSpectrum();
                        if (spectrum.getObservedMW() > Constants.maxPeptideMass) continue;

                        PGraph graph = spectrum.getPeakGraph();
                        spectrum.setCorrectedParentMW(graph.correctMW(dynamicPMCorrection));
                        TagPool tPool = (new TagSourceBuilder()).buildTagPool(spectrum);

                        DPHeap heatedPepts = (new OneMOD()).getHeatedPeptides(ixPDB, graph, tPool, considerIsotopeErr);
                        DPHeap tepidPepts = null;
                        if (Constants.maxPTMPerPeptide > 1) {
                            if (heatedPepts == null || !heatedPepts.isConfident()) {
                                tepidPepts = heatedPepts;
                                heatedPepts = (new MultiMOD()).getHeatedPeptides(ixPDB, graph, tPool, dynamicPMCorrection);
                            }
                        }

                        if (heatedPepts == null) continue;

                        final int numHeatedPeptides_f = numHeatedPeptides;
                        HeatedDB bitDB = (new DBHolder()).getHeatedDB(ixPDB, heatedPepts, tepidPepts,numHeatedPeptides_f);
                        TagTrie bitTrie = bitDB.getPartialDB(ixPDB);

                        tp = (new DynamicDBHolder()).dynamicMODeye(bitTrie, graph, tPool);
                        if (tp.size() > 0) {
                            if (candidates == null || candidates.get(0).compareTo(tp.get(0)) == 1) {
                                candidates = tp;
                                selected = JJ;
                            }
                        }
                    }
                    if (selected != -1) {
                        HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap = new HashMap<>();
                        try {
                            for (int k = 0; k < candidates.size(); k++) {
                                final int KK = k;
                                String tpSeq = candidates.get(KK).getPeptideSequence();
                                ArrayList<PeptideMatchToProtein> matchedProteins = seqToProtMap.get(tpSeq);

                                if (matchedProteins == null) {
                                    // 기존에 테이블에 없으면 새롭게 키 밸류 정보 등록 (최초)
                                    matchedProteins = ixPDB.getMatchProteins(tpSeq);
                                    seqToProtMap.put(tpSeq, matchedProteins);
                                }
                                // 기존에 있으면 테이블에 등록하지 않음. 즉 스킵합.
                            }
                        } finally {
                            results.put(scanIndex, new ResultEntry(scanTarget.get(selected), candidates, seqToProtMap));
                        }


                    }
                } catch (Exception e) {
                    //e.printStackTrace();
                } finally {
                    System.out.println(latch.getCount() + "<-" + totalSize);
                    latch.countDown();
                }
            });
        }


        latch.await();

        final String lastOutFileName = scaniter.getFileName();


        // Write results to file in original order
        try (PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(identifier + ".modplus.txt")))) {
            for (int i = 0; i < allScans.size(); i++) {
                ResultEntry entry = results.get(i);
                if (entry != null) {
                    out.println(">>" + lastOutFileName + "\t" + entry.scan.getHeader());

                    for (int k = 0; k < entry.candidates.size(); k++) {
                        AnsPeptide candidate = entry.candidates.get(k);
                        String tpSeq = candidate.getPeptideSequence();
                        ArrayList<PeptideMatchToProtein> matchedProteins = entry.seqToProtMap.get(tpSeq);
                        out.println(candidate.toMODPlus(entry.scan.getObservedMW(), matchedProteins));
                    }
                    out.println();
                }
            }
        }


        // ExecutorService 종료
        executor.shutdown();
        executor.awaitTermination(1, TimeUnit.HOURS);

        System.out.println("[MOD-Plus] Elapsed Time : " + (System.currentTimeMillis() - startTime) / 1000 + " Sec");
        return 0;
    }





    private static class ResultEntry {
        final MSMScan scan;
        final ArrayList<AnsPeptide> candidates;
        final HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap;

        ResultEntry(MSMScan scan, ArrayList<AnsPeptide> candidates,
                    HashMap<String, ArrayList<PeptideMatchToProtein>> seqToProtMap) {
            this.scan = scan;
            this.candidates = candidates;
            this.seqToProtMap = seqToProtMap;
        }
    }



}





