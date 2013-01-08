/*******************************************************************************
 * The ABAMS project
 * 
 * Copyright (c) 2012 University of British Columbia
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package ubic.BAMSandAllen;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices2;
import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.FocusedAnalysis.ExploreRegionNames;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.BAMSandAllen.geneFilters.PlaneRemoveFilter.Plane;
import ubic.basecode.dataStructure.StringToStringSetMap;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.util.FileTools;
import au.com.bytecode.opencsv.CSVReader;

public class RankedGeneListLoader {
    private static Log log = LogFactory.getLog( RankedGeneListLoader.class.getName() );
    static final int CORRELATION_TREND = 5;
    List<String> lines;
    String filename;
    Double threshold;

    public List<String> getLines() {
        return lines;
    }

    public void addToFront( List<String> addIn ) {
        lines.addAll( 0, addIn );
    }

    public RankedGeneListLoader( String filename ) throws Exception {
        this( filename, true );
    }

    public RankedGeneListLoader( String filename, boolean reverse ) throws Exception {
        this.filename = filename;
        lines = FileTools.getLines( filename );
        log.info( "Lines read:" + lines.size() + " " + filename );
        // chomp the last one if it's blank
        if ( lines.get( lines.size() - 1 ).equals( "" ) ) {
            lines.remove( lines.size() - 1 );
        }
        // assume the last is the first ranked gene
        if ( reverse ) reverse();
        threshold = null;
    }

    public RankedGeneListLoader( List<String> lines, String filename ) throws Exception {
        this.filename = filename;
        this.lines = lines;
        threshold = null;
    }

    /**
     * Writes information by using the matrix pair - a list of top genes and a the correlation after their removal
     * 
     * @param pair
     * @return
     * @throws Exception
     */
    public String writeCorrelationBasedData( MatrixPair pair ) throws Exception {
        return writeCorrelationBasedData( pair, true );
    }

    public String getResultFilename() {
        String result = filename.replace( "LOOGenesInOrder", "LOOResults" );
        return result;
    }

    public String getCorrelationLevelPoint() throws Exception {
        return writeEndSetGenes( true );
    }

    public String writeEndSetGenes( boolean keepSign ) throws Exception {
        List<String> topGenes = new LinkedList<String>();
        final int correlation_column = 1;
        final int gene_column = 2;
        String resultFilename = getResultFilename();
        CSVReader reader = new CSVReader( new FileReader( resultFilename ) );
        List<String[]> correlationLines = reader.readAll();
        double firstBaseLine = Double.parseDouble( correlationLines.get( 0 )[correlation_column] );

        boolean increase = firstBaseLine > 0;
        if ( !keepSign ) increase = !increase;

        List<String> rowsInFile = new LinkedList<String>();

        double apexCorrelation;

        if ( increase )
            apexCorrelation = Double.MIN_VALUE;
        else
            apexCorrelation = Double.MAX_VALUE;

        // find the apex
        for ( String[] line : correlationLines ) {
            double correlation = Double.parseDouble( line[correlation_column] );
            if ( increase && ( apexCorrelation < correlation ) ) {
                apexCorrelation = correlation;
            }
            if ( !increase && ( apexCorrelation > correlation ) ) apexCorrelation = correlation;
        }
        log.info( "Peak correlation:" + apexCorrelation );

        int count = 0;
        boolean reachedApex = false;
        // get genes after the apex, could use index instead
        for ( String[] line : correlationLines ) {
            count++;
            String row = line[gene_column];
            rowsInFile.add( row );
            double correlation = Double.parseDouble( line[correlation_column] );

            if ( correlation == apexCorrelation ) reachedApex = true;

            if ( reachedApex ) {
                topGenes.add( row );
            }
        }

        // the results file maybe missing genes at the end, add them in
        List<String> addToEnd = new LinkedList<String>( lines );
        addToEnd.removeAll( rowsInFile );
        log.info( "adding " + addToEnd.size() + " genes to end" );
        topGenes.addAll( addToEnd );

        threshold = ( ( double ) topGenes.size() ) / lines.size();
        // round to 5 decimals
        threshold = Double.parseDouble( String.format( "%.5g%n", threshold ) );

        String outFileName = filename + "." + topGenes.size() + "." + threshold + ".topGenes.txt";
        FileTools.stringsToFile( topGenes, new File( outFileName ) );
        log.info( "Wrote file: " + outFileName );

        return outFileName;
    }

    public void convertToNCBIIDs() throws Exception {
        ImageSeriesInfoLoader infoLoader = new ImageSeriesInfoLoader();
        List<String> IDs = new LinkedList<String>();
        int nullNCBI = 0;
        for ( String geneRowName : lines ) {
            Integer ncbiID = infoLoader.getNCBIIDFromRowName( geneRowName );
            if ( ncbiID != null ) {
                IDs.add( "" + ncbiID );
            } else {
                nullNCBI++;
            }
        }
        log.info( "Null NCBI genes:" + nullNCBI );
        lines = IDs;
    }

    @Deprecated
    public String writeCorrelationBasedData( MatrixPair pair, boolean keepSign ) throws Exception {
        // compute correlations
        // create pair? erg, no arguments
        // get correlation numbers
        reverse();

        double firstBaseLine = pair.getCorrelation( true );
        // make it more negative if it starts below zero
        boolean increase = firstBaseLine > 0;

        // if we are going against the current correlation sign - eg go from positive correlation to negative
        if ( !keepSign ) {
            increase = !increase;
        }

        // String result;// = "NumberRemoved\tCorrelation\n"+

        String newLine = System.getProperty( "line.separator" );
        String result = "";
        List<String> topGenes = new LinkedList<String>();

        // only work with the lines we have
        pair.setMatrixBDataRows( lines );
        long startTime = System.currentTimeMillis();

        result += "0\t" + pair.getCorrelation( true ) + "\t" + "AllGenes" + newLine;
        double lastCorrelation = pair.getCorrelation( true );
        int count = 0;
        for ( String row : lines ) {
            count++;
            pair.removeMatrixBDataRowFast( row );
            double correlation = pair.getCorrelation( true );
            result += count + "\t" + correlation + "\t" + row + newLine;
            if ( count % 1000 == 0 )
                log.info( "Removing:" + row + " Correlation after:" + correlation + " Removed:" + count );

            FileTools.stringToFile( count + "," + correlation + "," + row + "\n", new File( SetupParameters
                    .getDataFolder()
                    + "LOOResults." + startTime + ".txt" ), true );

            // if we are no longer gaining correlation then start keeping track
            // boolean decreasingCurrent = lastCorrelation - correlation > 0;

            // if we have more than four already, consider it a trend
            // if ( ( decreasingCurrent && increase ) || ( !decreasingCurrent && !increase )
            // || topGenes.size() > CORRELATION_TREND ) {
            // // if the non expressing genes were added to the end of the list then don't add this gene
            // // if ( count > 20000 ) {
            // topGenes.add( row );
            // log.info( "Top gene add " + count );
            // // }
            // } else {
            // // require a decrease for four genes in a row to start topGene list
            // if ( topGenes.size() != 0 && topGenes.size() < CORRELATION_TREND ) {
            // topGenes.clear();
            // log.info( "Top genes clear" );
            // }
            // }
            // lastCorrelation = correlation;
        }

        // // put it back in order
        // reverse();
        //
        // // make it so top gene is top
        // Collections.reverse( topGenes );
        //
        // FileTools.stringToFile( result, new File( filename + ".forPlotting" ) );
        // log.info( "Wrote file: " + filename + ".forPlotting" );
        //
        // FileTools.stringsToFile( topGenes, new File( filename + ".topGenes.txt" ) );
        // log.info( "Wrote file: " + filename + "." + topGenes.size() + ".topGenes.txt" );

        // RankedGeneListLoader tempLoader = new RankedGeneListLoader( topGenes, filename + "." + topGenes.size()
        // + ".topGenes" );
        // tempLoader.forErmineJ();

        return result;
    }

    public List<String> addMissing( boolean top ) throws Exception {
        List<String> originalRows = getOriginalRows();

        List<String> missing = new LinkedList<String>( originalRows );
        missing.removeAll( lines );
        // consistnent shuffling for now
        Collections.shuffle( missing, new Random( 1 ) );
        log.info( "Missing size:" + missing.size() );
        // log.info( "Missing:" + missing );
        // add to start of list
        for ( String missed : missing ) {
            if ( top )
                lines.add( 0, missed );
            else
                lines.add( missed );
        }
        return missing;
    }

    public List<String> addMissingToBottom() throws Exception {
        return addMissing( false );
    }

    public RankedGeneListLoader getShiftedList( RankedGeneListLoader list2 ) throws Exception {
        Map<String, Integer> shifts = new HashMap<String, Integer>();
        for ( int i = 0; i < lines.size(); i++ ) {
            String gene = lines.get( i );
            int j = list2.lines.indexOf( gene );
            int shift = i - j;
            shifts.put( gene, shift );
        }

        List<Map.Entry<String, Integer>> list = new LinkedList( shifts.entrySet() );
        Collections.sort( list, new Comparator() {
            public int compare( Object o1, Object o2 ) {
                return ( ( Comparable ) ( ( Map.Entry ) ( o1 ) ).getValue() ).compareTo( ( ( Map.Entry ) ( o2 ) )
                        .getValue() );
            }
        } );

        List<String> result = new LinkedList<String>();
        for ( Map.Entry<String, Integer> entry : list ) {
            // log.info( entry.getKey() + "->" + entry.getValue() );
            result.add( entry.getKey() );
        }
        return new RankedGeneListLoader( result, filename + ".shifted" );
    }

    public int getRank( String gene ) {
        return lines.indexOf( gene );
    }

    public void reverse() {
        Collections.reverse( lines );
    }

    public void removeRikStars() {
        List<String> result = new LinkedList<String>();
        int stars = 0;
        for ( String line : lines ) {
            if ( line.matches( ".*Rik[*].*" ) ) {
                line = line.replaceAll( "Rik[*]", "Rik" );
                stars++;
            }
            result.add( line );
        }
        log.info( "Stars removed:" + stars );
        lines = result;
    }

    public void removeImageIDs() {
        List<String> result = new LinkedList<String>();
        for ( String line : lines ) {
            line = line.replaceAll( "\\[.*\\]", "" );
            result.add( line );
        }
        lines = result;
    }

    public String getTabSepRanks() {
        String result = "gene\tscore\n";
        double rank = 0;
        for ( String line : lines ) {
            double scale = rank / lines.size();
            result += line + "\t" + scale + System.getProperty( "line.separator" );
            rank++;
        }
        return result;
    }

    public List<String> getRows() {
        return lines;
    }

    public int size() {
        return lines.size();
    }

    public String forErmineJ() throws Exception {
        removeRikStars();
        removeImageIDs();
        return writeOut();
    }

    public String forErmineJNCBI() throws Exception {
        convertToNCBIIDs();
        return writeOut( "NCBI" );
    }

    public String writeOut() throws Exception {
        return writeOut( "symbols" );
    }

    public String writeOut( String endfix ) throws Exception {
        log.info( "Size for writing:" + size() );
        String out_filename = filename + "." + endfix + ".forErmine";
        FileTools.stringToFile( getTabSepRanks(), new File( out_filename ) );
        log.info( "Wrote file: " + out_filename );
        return out_filename;
    }

    // public void makeChartData( MatrixPair pair ) throws Exception {
    // FileTools.stringToFile( getTabSepChartData( pair ), new File( filename + ".forPlotting" ) );
    // log.info( "Wrote file: " + filename + ".forPlotting" );
    // }

    public static List<String> getOriginalRows() throws Exception {
        AllenCatalogMatrices2 allenMatrices = new AllenCatalogMatrices2();
        DoubleMatrix<String, String> dataMatrix = allenMatrices.getFromDisk( "NewEnergies" );
        return dataMatrix.getRowNames();
    }

    public static void processDirectoryNoMissing( String directory, ConnectivityAndAllenExpressionMatrixPair pair )
            throws Exception {
        processDirectoryNoMissing( directory, pair, true );
    }

    public static void processDirectoryNoMissing( String directory, ConnectivityAndAllenExpressionMatrixPair pair,
            boolean keepSign ) throws Exception {
        Collection<File> files = FileTools.listDirectoryFiles( new File( directory ) );

        File genSetGenesFile = new File( directory + "/GeneSetVennMasterGenes.txt" );
        genSetGenesFile.delete();

        File genSetFile = new File( directory + "/GeneSetVennMaster.txt" );
        genSetFile.delete();

        for ( File file : files ) {
            String name = file.getName();
            if ( name.contains( "LOOGenesInOrder" ) && !name.contains( "forErmine" ) && !name.contains( "topGenes" )
                    && !name.contains( ".xls" ) ) {
                RankedGeneListLoader loaderIn = new RankedGeneListLoader( file.toString() );
                String topGenesFileName = loaderIn.writeEndSetGenes( keepSign );
                String SymbolsFilename = loaderIn.forErmineJ();

                // reload
                RankedGeneListLoader NCBIloaderIn = new RankedGeneListLoader( file.toString() );
                String NBCIFileName = NCBIloaderIn.forErmineJNCBI();

                // runErmineJORA( loaderIn.threshold, SymbolsFilename, false );
                runErmineJORA( loaderIn.threshold, NBCIFileName, true );

                // generate hascoronal list
                RankedGeneListLoader hasCoronal = new RankedGeneListLoader( file.toString() );
                hasCoronal.filename += ".hasCoronal";
                hasCoronal.retainHasCoronal();
                hasCoronal.forErmineJNCBI();

                // generate coronal list
                RankedGeneListLoader coronal = new RankedGeneListLoader( file.toString() );
                coronal.filename += ".coronal";
                coronal.retainCoronal();
                coronal.forErmineJNCBI();

                RankedGeneListLoader loaderTop = new RankedGeneListLoader( topGenesFileName );
                TopTenInfo topTen = new TopTenInfo( loaderTop, pair );
                topTen.writeExpressionInfo();

                String setName = name.substring( 16, name.indexOf( ".txt" ) );
                Util.addToVennMasterFile( genSetFile, loaderTop.lines, setName );

                // without IDS
                loaderTop.removeImageIDs();
                Util.addToVennMasterFile( genSetGenesFile, loaderTop.lines, setName );

            }
        }
    }

    public static void runErmineJORA( double threshold, String list_filename, boolean NCBI ) throws Exception {
        runErmineJORA( threshold, list_filename, NCBI, 5, 200 );
    }

    public static void runErmineJORA( double threshold, String list_filename, boolean NCBI, int minSize, int maxSize )
            throws Exception {
        // erg, so annoying
        List<String> commandParts = new LinkedList<String>();
        commandParts.add( SetupParameters.config.getString( "abams.ermineJ.bin" ) );
        if ( NCBI ) {
            commandParts.add( "--annots" );
            commandParts.add( "\"" + SetupParameters.config.getString( "abams.ermineJ.NCBIannots" ) + "\"" );
        } else {
            commandParts.add( "--annots" );
            commandParts.add( "\"" + SetupParameters.config.getString( "abams.ermineJ.annots" ) + "\"" );
        }
        commandParts.add( "--config" ); // config from GUI
        commandParts.add( "\"" + SetupParameters.config.getString( "abams.ermineJ.baseConfig" ) + "\"" );

        commandParts.add( "--classFile" );
        commandParts.add( "\"" + SetupParameters.config.getString( "abams.ermineJ.classFile" ) + "\"" );
        commandParts.add( "--output" );
        commandParts.add( "\"" + list_filename + ".ORA.tsv\"" );
        commandParts.add( "--threshold" );
        commandParts.add( "" + threshold );
        commandParts.add( "--test" );
        commandParts.add( "0" );
        commandParts.add( "--scoreFile" );
        commandParts.add( "\"" + list_filename + "\"" );
        commandParts.add( "--maxClassSize" );
        commandParts.add( ""  + maxSize);
        commandParts.add( "--minClassSize" );
        commandParts.add( "" + minSize );
        commandParts.add( "-j" ); // show genes

        // log.info( command );
        // execute!
        log.info( commandParts );
        Process p = Runtime.getRuntime().exec( ( String[] ) commandParts.toArray( new String[commandParts.size()] ) );
        log.info( "Waiting for ermineJ" );

        StreamGobbler outputGobbler = new StreamGobbler( p.getInputStream(), "ErmineJ" );
        outputGobbler.start();
        StreamGobbler errGobbler = new StreamGobbler( p.getErrorStream(), "ErmineJ" );
        errGobbler.start();
        p.waitFor();

    }

    static class StreamGobbler extends Thread {
        InputStream is;
        String type;

        StreamGobbler( InputStream is, String type ) {
            this.is = is;
            this.type = type;
        }

        public void run() {
            try {
                InputStreamReader isr = new InputStreamReader( is );
                BufferedReader br = new BufferedReader( isr );
                String line = null;
                while ( ( line = br.readLine() ) != null )
                    System.out.println( type + ">" + line );
            } catch ( IOException ioe ) {
                ioe.printStackTrace();
            }
        }
    }

    @Deprecated
    public static void processDirectory( String directory, boolean hasNonExp ) throws Exception {
        Collection<File> files = FileTools.listDirectoryFiles( new File( directory ) );
        for ( File file : files ) {
            RankedGeneListLoader loaderIn = new FullRankedGeneListLoader( file.toString(), hasNonExp );
            loaderIn.forErmineJ();
        }
    }

    /**
     * @param args
     */

    public String getTotalExpresion() throws Exception {
        String result = "";
        AllenCatalogMatrices2 allenMatrices = new AllenCatalogMatrices2();
        DoubleMatrix<String, String> dataMatrix = allenMatrices.getFromDisk( "NewEnergies" );

        for ( String geneRow : lines ) {
            double[] values = dataMatrix.getRowByName( geneRow );
            double sum = 0;
            for ( double d : values )
                sum += d;
            result += geneRow + "\t" + sum + System.getProperty( "line.separator" );
        }

        return result;
    }

    public void retainHasCoronal() throws Exception {
        ImageSeriesInfoLoader loader = new ImageSeriesInfoLoader();
        List<String> newLines = new LinkedList<String>();
        for ( String row : lines ) {
            boolean hasCoronal = loader.hasCoronalImageFromRowName( row );
            if ( hasCoronal ) newLines.add( row );
        }
        this.lines = newLines;
    }

    public void retainCoronal() throws Exception {
        ImageSeriesInfoLoader loader = new ImageSeriesInfoLoader();
        List<String> newLines = new LinkedList<String>();
        for ( String row : lines ) {
            String planeString = loader.getPlaneFromRowName( row );
            if ( planeString.equals( "coronal" ) ) newLines.add( row );
        }
        this.lines = newLines;
    }

    public static void main( String[] args ) throws Exception {

        // RankedGeneListLoader jesseIdea = new RankedGeneListLoader(
        // "/grp/java/workspace/BAMSandAllen/data/rankedGenes/ranked once/Incoming.partialcon.genes.txt" );

        boolean removeNonExp = true;
        boolean useVirtual = true;
        boolean keepSign = true;

        // ERG, gene info top ten results varies for direction
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;

        ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityAndExpression(
                direction, useVirtual, removeNonExp );
        // ConnectivityAndAllenExpressionMatrixPair pair = ExpressionMatrixPairFactory.connectivityPartial( direction,
        // false, RegressMatrix.CONNECTIVITY, useVirtual, removeNonExp, true );

        // need a partialcon pair
        // jesseIdea.writeCorrelationBasedData( pair );
        // System.exit( 1 );
        // pair.applyGeneFilter( new PrefixGeneFilter( "Drd" ) );
        // log.info( "Genes:" + Util.getUniqueGenes( pair.getMatrixBDataRows() ) );
        // log.info( pair.getCorrelation() );
        // pair.writeImages();
        //
        // System.exit( 1 );

        // processDirectoryNoMissing( "/grp/java/workspace/BAMSandAllen/data/rankedGenes/direct increasing", pair );

        // processDirectoryNoMissing( "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final nobed", pair );
        processDirectoryNoMissing( "/grp/java/workspace/BAMSandAllen/data/rankedGenes/LitCon", pair );
        System.exit( 1 );

        // LOOGenesInOrder.out.partialcon.txt.329.0.014448.topGenes.txt
        // LOOGenesInOrder.in.partialcon.txt.410.0.018005.topGenes.txt
        // LOOGenesInOrder.space.7.txt.435.0.019093.topGenes.txt

        /*
         * RankedGeneListLoader aLook = new RankedGeneListLoader(
         * "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final
         * nobed/LOOGenesInOrder.space.7.txt.435.0.019093.topGenes.txt" );
         */
        // ConnectivityAndAllenExpressionMatrixPair.NewEnergies.incoming.rOrder
        // ConnectivityAndAllenExpressionMatrixPair.NewEnergies.space.rOrder
        // ConnectivityAndAllenExpressionMatrixPair.NewEnergies.outgoing.rOrder
        // LOOGenesInOrder.out.partialcon.txt.329.0.014448.topGenes.txt
        // LOOGenesInOrder.in.partialcon.txt.410.0.018005.topGenes.txt
        // LOOGenesInOrder.space.7.txt.435.0.019093.topGenes.txt
        RankedGeneListLoader aLook = new RankedGeneListLoader(
                "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final nobed/LOOGenesInOrder.out.partialcon.txt.329.0.014448.topGenes.txt" );
        pair.setMatrixBDataRows( aLook.lines );
        // log.info( pair.getCorrelation() );
        // pair.writeRMatrices();
        // pair.writeImages();
        pair.orderDataRows( aLook.lines );
        log.info( pair.getCorrelation() );
        System.exit( 1 );
        pair.writeRMatrices();
        pair.writeImages();

        ExploreRegionNames explore = new ExploreRegionNames( pair );
        StringToStringSetMap parents = explore.getParents();

        String focusRegion = "Midbrain";
        Set<String> ROIs = parents.get( focusRegion );
        // some may have no exp
        ROIs.retainAll( pair.getAllenDataColNames() );
        pair.removeAllenCols( ROIs );
        // Hindbrain
        // Interbrain
        // Midbrain
        // Cerebrum
        pair.run();

        log.info( focusRegion + ":" + pair.getCorrelation() );
        pair.test( 1000 );

        System.exit( 1 );

        List<String> lines = aLook.lines;
        lines.remove( 0 );
        lines.remove( 1 );
        lines.remove( 2 );
        lines.remove( 3 );
        lines.remove( 4 );
        lines.remove( 5 );
        lines.remove( 6 );
        lines.remove( 7 );
        lines.remove( 8 );
        lines.remove( 9 );
        pair.setMatrixBDataRows( aLook.lines );
        log.info( pair.getCorrelation() );
        log.info( lines.size() );

        // processDirectoryNoMissing( "/grp/java/workspace/BAMSandAllen/data/rankedGenes/for Ray", pair );
        // processDirectoryNoMissing( "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final 962 regions", pair
        // );
        // processDirectoryNoMissing( "/grp/java/workspace/BAMSandAllen/data/rankedGenes/decreasing after lab meeting",
        // pair, keepSign );

    }
}
