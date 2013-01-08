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
package ubic.BAMSandAllen.FocusedAnalysis;

import static ubic.BAMSandAllen.Util.intersect;
import static ubic.BAMSandAllen.Util.intersectSize;
import static ubic.BAMSandAllen.Util.subtract;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes;
import ubic.BAMSandAllen.AllenDataLoaders.AllenAtlasAnnotationLoader;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.AnalyzeBAMSandAllenGenes.Direction;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.MatrixPairs.ConnectivityAndAllenExpressionMatrixPair;
import ubic.BAMSandAllen.MatrixPairs.MatrixPair;
import ubic.basecode.dataStructure.StringToStringSetMap;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class ExploreRegionNames {
    private static Log log = LogFactory.getLog( ExploreRegionNames.class.getName() );

    MatrixPair pair;
    List<String> aNames;
    List<String> bNames;
    StructureCatalogLoader loader;

    public ExploreRegionNames( MatrixPair pair ) throws Exception {
        this.pair = pair;
        aNames = pair.getMatrixA().getColNames();
        bNames = new LinkedList<String>( pair.getMatrixB().getColNames() );
        log.info( bNames.size() );
        log.info( bNames );
        log.info( "Removing bed nuclei manually" );
        bNames.remove( "Bed nuclei of the stria terminalis" );
        loader = new StructureCatalogLoader();
    }

    public void printSizes( Collection<String> allenRegions ) throws Exception {
        AllenAtlasAnnotationLoader spaceLoader = new AllenAtlasAnnotationLoader();
        DoubleMatrix<String, String> volumes = spaceLoader.getVolumeMatrix();
        for ( String region : allenRegions ) {
            System.out.println( region + "," + volumes.getByKeys( "volume", region ) );
        }
    }

    public StringToStringSetMap getParents() {
        StringToStringSetMap parents = new StringToStringSetMap();
        for ( String allenRegion : loader.getRegions() ) {
            Set<String> regionParents = loader.getParents( allenRegion );
            for ( String parent : regionParents ) {
                parents.put( parent, allenRegion );
            }
        }
        return parents;
    }

    public Collection<String> loader() throws Exception {
        StringToStringSetMap parents = getParents();

        Set<String> nullBAMS = new HashSet<String>();
        Set<String> notLeaf = new HashSet<String>();
        Set<String> leafs = new HashSet<String>();
        Set<String> nucleus = new HashSet<String>();
        Set<String> goodMapping = new HashSet<String>();
        Set<String> hasConnectivityData = new HashSet<String>();
        for ( String allenRegion : loader.getRegions() ) {
            Set<String> aMappedNames = pair.convertBNametoA( allenRegion );
            if ( aMappedNames != null && !aMappedNames.isEmpty() ) {
                if ( intersectSize( aNames, aMappedNames ) > 0 ) {
                    hasConnectivityData.add( allenRegion );
                }
            }

            if ( !loader.getLeafs().contains( allenRegion ) ) {
                notLeaf.add( allenRegion );
            } else {
                leafs.add( allenRegion );
            }
            Set<String> result = loader.getBAMSMappedRegions( allenRegion );
            if ( result == null ) {
                nullBAMS.add( allenRegion );
            } else {
                goodMapping.add( allenRegion );
            }
            if ( allenRegion.toLowerCase().contains( "nucleus" ) ) {
                nucleus.add( allenRegion );
            }

        }

        log.info( "nucleus:" + nucleus.size() );
        log.info( "leafs:" + leafs.size() );
        log.info( "notLeaf:" + notLeaf.size() );
        log.info( "Cortex is not a leaf?" + notLeaf.contains( "Cerebral cortex" ) );
        log.info( "Cortex children:" + loader.getChildren( "Cerebral cortex" ) );
        log.info( "Ammons horn is not a leaf?" + notLeaf.contains( "Ammon\'s Horn" ) );
        log.info( "Ammons children:" + loader.getChildren( "Ammon\'s Horn" ) );
        log
                .info( "Posterodorsal preoptic nucleus is not a leaf?"
                        + notLeaf.contains( "Posterodorsal preoptic nucleus" ) );
        log.info( "Posterodorsal preoptic nucleus children?" + notLeaf.contains( "Posterodorsal preoptic nucleus" ) );

        log.info( "Cortical plate children:" + loader.getChildren( "Cortical plate" ) );
        log.info( "nullBAMS:" + nullBAMS.size() );
        log.info( "goodMapping:" + goodMapping.size() );
        log.info( "In loader but have no expression:" + subtract( loader.getRegions(), bNames ) );

        log.info( "Pair amount (those with data):" + bNames.size() );
        log.info( "Loader regions:" + loader.getRegions().size() );
        log.info( loader.getRegions().contains( "Nucleus ambiguus" ) );
        log.info( "Regions missing matrixB cols(exp data):" + subtract( loader.getRegions(), bNames ).size() );
        log.info( "Regions missing matrixB cols(exp data):" + subtract( loader.getRegions(), bNames ) );
        log.info( "Leafs with no mapping:" + intersectSize( leafs, nullBAMS ) );
        log.info( "Leafs with no mapping:" + intersect( leafs, nullBAMS ) );
        log.info( "Mapped leafes with data:" + intersectSize( leafs, goodMapping, bNames ) );

        Collection<String> good = ( Collection<String> ) intersect( leafs, hasConnectivityData, bNames );
        log.info( "in good:" + good.contains( "Nucleus ambiguus" ) );
        log.info( "Final set:" + intersectSize( hasConnectivityData, good ) );
        log.info( "Final set minus leafs with data:"
                + ( intersectSize( hasConnectivityData, good ) - intersectSize( leafs, goodMapping, bNames ) ) );
        log.info( "nucleus and leafes with data:" + intersectSize( leafs, nucleus, bNames ) );
        log.info( "using cols that are nucleus:" + intersectSize( good, nucleus ) );
        log.info( "kept regions that are not nucleus:" + subtract( good, nucleus ) );

        // log.info( "Midbrain:" + parents.get( "Midbrain" ).size() );
        // log.info( "Hindbrain:" + parents.get( "Hindbrain" ).size() );
        // log.info( "Interbrain:" + parents.get( "Interbrain" ).size() );
        // log.info( "Cerebrum:" + parents.get( "Cerebrum" ).size() );
        // log.info( "Cerebellum:" + parents.get( "Cerebellum" ).size() );
        log.info( "nucleus in final:" + intersectSize( nucleus, good ) );

        log.info( "Midbrain kept:" + intersectSize( parents.get( "Midbrain" ), good ) + " before:"
                + parents.get( "Midbrain" ).size() );
        log.info( "Hindbrain kept:" + intersectSize( parents.get( "Hindbrain" ), good ) + " before:"
                + parents.get( "Hindbrain" ).size() );
        log.info( "Interbrain kept:" + intersectSize( parents.get( "Interbrain" ), good ) + " before:"
                + parents.get( "Interbrain" ).size() );
        log.info( "Cerebrum kept:" + intersectSize( parents.get( "Cerebrum" ), good ) + " before:"
                + parents.get( "Cerebrum" ).size() );
        log.info( "Cerebellum kept:" + intersectSize( parents.get( "Cerebellum" ), good ) + " before:"
                + parents.get( "Cerebellum" ).size() );
        log.info( "Cerebellum:" + intersect( good, parents.get( "Cerebellum" ) ) );

        AllenAtlasAnnotationLoader atlas = new AllenAtlasAnnotationLoader();
        log.info( "Total Volume all atlas regions:" + atlas.getTotalVolume( null ) );
        log.info( "Total Volume all atlas exp:" + atlas.getTotalVolume( bNames ) );
        log.info( "Volume contained in nucleus regions:" + atlas.getTotalVolume( nucleus ) );
        log.info( "Volume contained in final regions:" + atlas.getTotalVolume( good ) );
        //printSizes( good );
        return good;
    }

    public void printRelations() {
        log.info( "Bnames size:" + bNames.size() );
        log.info( "Anames size:" + aNames.size() );
        log.info( bNames.contains( "Nucleus y" ) );
        log.info( aNames.contains( "Nucleus y" ) );
        for ( String matrixAColumn : aNames ) {

            Set<String> matrixBColumns = pair.convertANametoB( matrixAColumn );
            // we may still have been mapped to a brain region with no data
            // requires mapping and matrix data!
            int BSizeBefore = matrixBColumns.size();
            Set<String> removed = new HashSet<String>( matrixBColumns );
            removed.removeAll( bNames );
            matrixBColumns.retainAll( bNames );

            int bSizeDiff = BSizeBefore - matrixBColumns.size();
            if ( bSizeDiff != 0 ) {
                // log.info( "B size changed by (no connect info):" + bSizeDiff );
                log.info( "No connection info for:   " + matrixAColumn + "->" + removed );
            }

            if ( matrixBColumns.size() > 1 ) {
                log.info( "Merging " + matrixBColumns + " into " + matrixAColumn );
            }
        }

    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        Direction direction = AnalyzeBAMSandAllenGenes.Direction.ANYDIRECTION;
        boolean squareMatrix = false;
        ConnectivityAndAllenExpressionMatrixPair forR = new ConnectivityAndAllenExpressionMatrixPair(
                new BrainRegionClassSelector(), true, false, squareMatrix, Double.NaN, "NewEnergies", direction );
        ExploreRegionNames test = new ExploreRegionNames( forR );
        test.printRelations();
        //forR.printMappingRelations();
        test.loader();
    }

}
