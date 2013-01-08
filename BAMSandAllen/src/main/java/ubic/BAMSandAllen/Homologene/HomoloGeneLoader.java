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
package ubic.BAMSandAllen.Homologene;

import java.io.FileReader;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.ABAMSDataMatrix;
import ubic.BAMSandAllen.RankedGeneListLoader;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.ImageSeriesInfoLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.FocusedAnalysis.CellTypes;
import ubic.BAMSandAllen.MatrixPairs.ExpressionMatrixPairFactory;
import ubic.BAMSandAllen.geneFilters.GeneFilter;
import ubic.BAMSandAllen.geneFilters.PlaneRemoveFilter;
import ubic.basecode.dataStructure.CountingMap;
import ubic.basecode.dataStructure.StringToStringSetMap;
import ubic.basecode.dataStructure.params.ParamKeeper;
import au.com.bytecode.opencsv.CSVReader;

//Given a NCBI gene for mouse gene
//
//map from mouse NCBI to other species
//
//has homolog in Saccharomyces cerevisiae?
//has homolog gene is in c.elegans?
//has homolog gene is in fly?
// 
// should be refactored to use basedCode/Gemma homologene class

public class HomoloGeneLoader {

    private static Log log = LogFactory.getLog( CellTypes.class.getName() );

    StringToStringSetMap clusterIDtoSpecies;
    Map<String, String> mouseIDtoCluster, humanIDtoCluster, elegansSymboltoCluster, humanSymboltoCluster;
    final int CLUSTERID = 0; // HomoloGene group id
    final int SPECIESID = 1;
    final int GENEID = 2;
    final int SYMBOL = 3;

    final String MOUSESPECIES = "10090";
    final String HUMANSPECIES = "9606";
    final String YEASTSPECIES = "4932";
    final String FLYSPECIES = "7227";
    final String WORMSPECIES = "6239";
    final String RATSPECIES = "10116";

    public HomoloGeneLoader() throws Exception {
        clusterIDtoSpecies = new StringToStringSetMap();
        mouseIDtoCluster = new HashMap<String, String>();
        elegansSymboltoCluster = new HashMap<String, String>();
        humanSymboltoCluster = new HashMap<String, String>();
        humanIDtoCluster = new HashMap<String, String>();

        String filename = SetupParameters.config.getString( "abams.homologene.database" );

        FileReader fin = new FileReader( filename );
        CSVReader reader = new CSVReader( fin, '\t' );
        String[] line;

        while ( ( line = reader.readNext() ) != null ) {
            String cluster = line[CLUSTERID];
            String species = line[SPECIESID];
            String geneid = line[GENEID];
            String symbol = line[SYMBOL];

            clusterIDtoSpecies.put( cluster, species );

            if ( species.equals( MOUSESPECIES ) ) {
                mouseIDtoCluster.put( geneid, cluster );
            }
            if ( species.equals( WORMSPECIES ) ) {
                elegansSymboltoCluster.put( symbol, cluster );
            }
            if ( species.equals( HUMANSPECIES ) ) {
                humanSymboltoCluster.put( symbol, cluster );
                humanIDtoCluster.put( geneid, cluster );
            }
        }
        fin.close();

    }

    public Set<String> getSpecies( String ClusterID ) {
        return clusterIDtoSpecies.get( ClusterID );
    }

    public String getMouseIDFromClusterID( String clusterID ) {
        for ( String key : mouseIDtoCluster.keySet() ) {
            if ( mouseIDtoCluster.get( key ).equals( clusterID ) ) return key;
        }
        return null;
    }

    public String getClusterIDFromWormSymbol( String sym ) {
        return elegansSymboltoCluster.get( sym );
    }

    public Set<String> getHumanSymbolFromMouseID( Collection<String> mouse ) throws Exception {
        Set<String> result = new HashSet<String>();
        for ( String mouseID : mouse ) {
            result.addAll( getHumanSymbolFromMouseID( mouseID ) );
        }
        return result;
    }

    public Set<String> getHumanSymbolFromMouseID( String mouse ) throws Exception {
        return getHumanIDorSymFromMouseID( mouse, humanSymboltoCluster );
    }

    public Set<String> getHumanIDsFromMouseID( Collection<String> mouse ) throws Exception {
        Set<String> result = new HashSet<String>();
        for ( String mouseID : mouse ) {
            result.addAll( getHumanIDsFromMouseID( mouseID ) );
        }
        return result;
    }

    public Set<String> getHumanIDsFromMouseID( String mouse ) throws Exception {
        return getHumanIDorSymFromMouseID( mouse, humanIDtoCluster );
    }

    public Set<String> getHumanIDorSymFromMouseID( String mouse, Map<String, String> map ) throws Exception {
        String clusterOfMouse = mouseIDtoCluster.get( mouse );
        Set<String> result = new HashSet<String>();
        for ( String humanIDorSym : map.keySet() ) {
            if ( map.get( humanIDorSym ).equals( clusterOfMouse ) ) result.add( humanIDorSym );
        }
        return result;
    }

    public Set<String> getSpeciesFromMouseID( String mouseID ) {
        return clusterIDtoSpecies.get( mouseIDtoCluster.get( mouseID ) );
    }

    public boolean hasHomologCluster( String mouseID ) {
        return getSpeciesFromMouseID( mouseID ) != null;
    }

    public boolean hasHomolog( String speciesID, String mouseID ) {
        if ( hasHomologCluster( mouseID ) ) {
            Set<String> species = getSpeciesFromMouseID( mouseID );
            return species.contains( speciesID );
        } else
            return false;
    }

    public boolean hasHumanHomolog( String mouseID ) {
        return hasHomolog( HUMANSPECIES, mouseID );
    }

    public boolean hasYeastHomolog( String mouseID ) {
        return hasHomolog( YEASTSPECIES, mouseID );
    }

    public boolean hasFlyHomolog( String mouseID ) {
        return hasHomolog( FLYSPECIES, mouseID );
    }

    public boolean hasWormHomolog( String mouseID ) {
        return hasHomolog( WORMSPECIES, mouseID );
    }

    public boolean hasRatHomolog( String mouseID ) {
        return hasHomolog( RATSPECIES, mouseID );
    }

    public CountingMap<String> getHomologCounts( Collection<String> allenRows ) throws Exception {
        ImageSeriesInfoLoader imageInfo = new ImageSeriesInfoLoader();
        Set<String> NCBIDS = new HashSet<String>();
        for ( String row : allenRows ) {
            String NCBID = imageInfo.getNCBIIDFromRowName( row ) + "";
            NCBIDS.add( NCBID );
        }
        return getHomologCountsFromNCBIDS( NCBIDS );
    }

    public CountingMap<String> getHomologCountsFromNCBIDS( Collection<String> NCBIDS ) throws Exception {
        CountingMap<String> result = new CountingMap<String>();
        for ( String NCBID : NCBIDS ) {
            result.increment( "totalCount" );
            if ( !hasHomologCluster( NCBID ) ) result.increment( "not in homologene" );
            if ( hasFlyHomolog( NCBID ) ) result.increment( "Fly" );
            if ( hasRatHomolog( NCBID ) ) result.increment( "Rat" );
            if ( hasWormHomolog( NCBID ) ) result.increment( "Worm" );
            if ( hasYeastHomolog( NCBID ) ) result.increment( "Yeast" );
            if ( hasHumanHomolog( NCBID ) ) result.increment( "Human" );
            if ( hasWormHomolog( NCBID ) || hasYeastHomolog( NCBID ) || hasFlyHomolog( NCBID ) )
                result.increment( "Worm or Yeast or Fly" );
        }
        return result;
        // geneStats.put( "NCBI ID", "" + imageInfo.getNCBIIDFromRowName( geneRow ) );

    }

    /**
     * @param args
     */
    public static Set<String> getKaufmanMantel() {
        Set<String> wormGenes = new HashSet<String>();
        wormGenes.add( "che-2" );
        wormGenes.add( "daf-28" );
        wormGenes.add( "egl-44" );
        wormGenes.add( "mec-12" );
        wormGenes.add( "mgl-2" );
        wormGenes.add( "mps-1" );
        wormGenes.add( "pef-1" );
        wormGenes.add( "slt-1" );
        wormGenes.add( "sup-9" );
        wormGenes.add( "unc-30" );
        wormGenes.add( "unc-5" );
        wormGenes.add( "unc-93" );
        wormGenes.add( "ace-1" );
        wormGenes.add( "ace-2" );
        wormGenes.add( "bra-1" );
        wormGenes.add( "cat-1" );
        wormGenes.add( "cat-2" );
        wormGenes.add( "ceh-23" );
        wormGenes.add( "che-3" );
        wormGenes.add( "daf-28" );
        wormGenes.add( "dat-1" );
        wormGenes.add( "deg-1" );
        wormGenes.add( "dop-2" );
        wormGenes.add( "egl-3" );
        wormGenes.add( "glr-8" );
        wormGenes.add( "gpa-14" );
        wormGenes.add( "gpa-16" );
        wormGenes.add( "gpa-2" );
        wormGenes.add( "gpa-3" );
        wormGenes.add( "gpa-9" );
        wormGenes.add( "hcp-3" );
        wormGenes.add( "ida-1" );
        wormGenes.add( "ina-1" );
        wormGenes.add( "kal-1" );
        wormGenes.add( "kin-29" );
        wormGenes.add( "kvs-1" );
        wormGenes.add( "lin-11" );
        wormGenes.add( "mab-9" );
        wormGenes.add( "mdl-1" );
        wormGenes.add( "mec-7" );
        wormGenes.add( "nlp-11" );
        wormGenes.add( "nlp-6" );
        wormGenes.add( "osm-3" );
        wormGenes.add( "osm-9" );
        wormGenes.add( "pef-1" );
        wormGenes.add( "rig-1" );
        wormGenes.add( "sem-4" );
        wormGenes.add( "ser-4" );
        wormGenes.add( "sra-10" );
        wormGenes.add( "sra-11" );
        wormGenes.add( "sre-37" );
        wormGenes.add( "sup-9" );
        wormGenes.add( "tax-2" );
        wormGenes.add( "tax-4" );
        wormGenes.add( "tol-1" );
        wormGenes.add( "tph-1" );
        wormGenes.add( "ttx-3" );
        wormGenes.add( "unc-30" );
        wormGenes.add( "unc-5" );
        wormGenes.add( "unc-53" );
        wormGenes.add( "unc-6" );
        wormGenes.add( "unc-8" );
        wormGenes.add( "unc-93" );
        wormGenes.add( "vab-15" );
        return wormGenes;
    }

    public static Set<String> getKaufmanBoth() {
        Set<String> wormGenes = new HashSet<String>();
        wormGenes.add( "che-2" );
        wormGenes.add( "mgl-2" );
        wormGenes.add( "mps-1" );
        wormGenes.add( "pef-1" );
        wormGenes.add( "unc-5" );
        wormGenes.add( "ceh-23" );
        wormGenes.add( "che-3" );
        wormGenes.add( "gpa-3" );
        wormGenes.add( "kin-29" );
        wormGenes.add( "kvs-1" );
        wormGenes.add( "lin-11" );
        wormGenes.add( "osm-3" );
        wormGenes.add( "osm-9" );
        wormGenes.add( "tax-2" );
        wormGenes.add( "tax-4" );
        return wormGenes;
    }

    public static void runCElegans() throws Exception {
        ImageSeriesInfoLoader geneInfo = new ImageSeriesInfoLoader();
        HomoloGeneLoader loader = new HomoloGeneLoader();

        Set<String> wormGenes = getKaufmanMantel();
        Set<String> mouseGenes = new HashSet<String>();

        for ( String wormGene : wormGenes ) {
            String clusterID = loader.getClusterIDFromWormSymbol( wormGene );
            if ( clusterID == null ) {
                log.info( wormGene + " null cluster" );
            } else {
                Set<String> species = loader.getSpecies( clusterID );
                if ( species.contains( loader.MOUSESPECIES ) ) {
                    String mouseID = loader.getMouseIDFromClusterID( clusterID );
                    String mouseSymbol = geneInfo.getSymbolFromNCBIID( mouseID );
                    log.info( wormGene + " hit, mouseID:" + mouseID + " " + mouseSymbol );
                    if ( mouseSymbol != null ) mouseGenes.add( mouseSymbol );
                } else {
                    log.info( wormGene + " miss" );
                }
            }
        }
        for ( String mouseGene : mouseGenes ) {
            System.out.println( mouseGene + "\tKaufmanMantel" );
        }
    }

    public static void runABAMSsets() throws Exception {
        // TODO Auto-generated method stub
        HomoloGeneLoader loader = new HomoloGeneLoader();
        log.info( loader.getSpeciesFromMouseID( "320405" ) );
        log.info( loader.hasHumanHomolog( "320405" ) );
        log.info( loader.hasWormHomolog( "22138" ) );
        log.info( loader.hasRatHomolog( "22138" ) );
        log.info( loader.hasFlyHomolog( "22138" ) );
        log.info( loader.hasYeastHomolog( "22138" ) );
        boolean removeNonExp = true;
        ABAMSDataMatrix matrix = ExpressionMatrixPairFactory.getEnergyMatrix( Direction.OUTGOING, removeNonExp );
        // log.info( loader.getHomologCounts( matrix.getRowNames() ) );

        ParamKeeper keeper = new ParamKeeper();
        Map<String, String> params;

        params = convertIntMap( loader.getHomologCounts( matrix.getRowNames() ) );
        params.put( "Name", "All" );
        keeper.addParamInstance( params );

        String filename = "/grp/java/workspace/BAMSandAllen/data/rankedGenes/near final ammon/";
        Set<String> files = new HashSet<String>();
        files.add( filename + "LOOGenesInOrder.in.partialcon.ammon.txt.452.0.01985.topGenes.txt" );
        files.add( filename + "LOOGenesInOrder.out.partialcon.ammon.txt.374.0.016424.topGenes.txt" );
        files.add( filename + "LOOGenesInOrder.space.ammon.txt.420.0.018435.topGenes.txt" );

        for ( String file : files ) {
            RankedGeneListLoader aLook = new RankedGeneListLoader( file );
            params = convertIntMap( loader.getHomologCounts( aLook.getRows() ) );
            params.put( "Name", file );
            keeper.addParamInstance( params );
        }

        // axon guidance
        //
        String file = "/grp/java/workspace/BAMSandAllen/data/G2Cdbgenelists/AxonGuidanceGOGroup.txt";
        params = loader.getCountsFromGeneNameFile( matrix, file );
        keeper.addParamInstance( params );

        file = "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/top100-pattern-a.txt";
        RankedGeneListLoader aLook = new RankedGeneListLoader( file );
        params = convertIntMap( loader.getHomologCounts( aLook.getRows() ) );
        params.put( "Name", file );
        keeper.addParamInstance( params );

        file = "/grp/java/workspace/BAMSandAllen/data/ABACellEnriched/patterns/top100-pattern-b.txt";
        aLook = new RankedGeneListLoader( file );
        params = convertIntMap( loader.getHomologCounts( aLook.getRows() ) );
        params.put( "Name", file );
        keeper.addParamInstance( params );

        // coronal
        GeneFilter filter = new PlaneRemoveFilter( PlaneRemoveFilter.Plane.SAGITTAL );
        List<String> removeRows = filter.getRowsToRemove( matrix );
        log.info( "Removed:" + removeRows.size() + " Name:" + filter.getName() );
        matrix = matrix.removeRows( removeRows );
        params = convertIntMap( loader.getHomologCounts( matrix.getRowNames() ) );
        params.put( "Name", "Coronal" );
        keeper.addParamInstance( params );

        keeper.writeExcel( SetupParameters.getDataFolder() + "homologene.xls" );
    }

    private Map<String, String> getCountsFromGeneNameFile( ABAMSDataMatrix matrix, String file ) throws Exception {
        Map<String, String> params;
        Set<String> NCBIDs = Util.getNCBIDSFromGeneNames( file, matrix.getRowNames() );
        params = convertIntMap( getHomologCountsFromNCBIDS( NCBIDs ) );
        params.put( "Name", file );
        return params;
    }

    // uggly!
    public static Map<String, String> convertIntMap( Map<String, Integer> in ) {
        Map<String, String> result = new HashMap<String, String>();
        for ( String key : in.keySet() ) {
            result.put( key, "" + in.get( key ) );
        }
        return result;
    }

    public static void main( String[] args ) throws Exception {
        HomoloGeneLoader loader = new HomoloGeneLoader();
        log.info( loader.getHumanSymbolFromMouseID( "320405" ).iterator().next() );

        // String endFile = "LOOGenesInOrder.outgoing.partialcon.txt";
        String endFile = "LOOGenesInOrder.outgoing.partialcon.txt.topGenes.txt";
        String baseFolder = SetupParameters.config.getString( "abams.dataFolder" ) + "rankedGenes/near final ammon/";
        // String outputFile = SetupParameters.config.getString( "abams.dataFolder" ) + "/topten/" + endFile;
        String outputFile = baseFolder + endFile;

        // rankedLoader.get
        // Util.getNCBIDSFromGeneNames

        // runABAMSsets();
    }
}
