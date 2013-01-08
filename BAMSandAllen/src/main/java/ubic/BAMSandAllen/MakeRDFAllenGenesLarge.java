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

import static ubic.BAMSandAllen.Vocabulary.getBAMSURI;
import static ubic.BAMSandAllen.Vocabulary.getBrainLinksURI;
import static ubic.BAMSandAllen.Vocabulary.hasAllen17Mapping;
import static ubic.BAMSandAllen.Vocabulary.hasAllenEnrichedGene;
import static ubic.BAMSandAllen.Vocabulary.hasStructureCatalogMapping;

import java.io.FileWriter;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.AllenDataLoaders.AllenCatalogMatrices;
import ubic.BAMSandAllen.AllenDataLoaders.AllenMajorMatrices;
import ubic.BAMSandAllen.AllenDataLoaders.AllenTop50DataLoader;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.BAMSDataLoaders.BAMSDataLoader;

import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.rdf.model.Model;
import com.hp.hpl.jena.rdf.model.ModelFactory;
import com.hp.hpl.jena.rdf.model.Resource;

public class MakeRDFAllenGenesLarge {
    protected static Log log = LogFactory.getLog( MakeRDFAllenGenesLarge.class );

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // load the mappings - CSV

        // mapping between BAMS and Allen
        Map<String, Set<String>> BAMStoAllenMapping = new FullMappingLoader().getAllentoBamsMapping();

        // BAMS
        BAMSDataLoader BAMSData = new BAMSDataLoader();

        // Allen
        AllenCatalogMatrices catalogExpression = new AllenCatalogMatrices();
        AllenMajorMatrices allenMatrices = new AllenMajorMatrices();
        List allen17Names = allenMatrices.getLevels().getColNames();

        AllenTop50DataLoader allenLoader = new AllenTop50DataLoader();

        StructureCatalogLoader allenStructureLoader = new StructureCatalogLoader();

        Model allenDataRDFModel = ModelFactory.createDefaultModel();
        allenDataRDFModel.setNsPrefix( "brianLink", getBrainLinksURI() );
        allenDataRDFModel.setNsPrefix( "BAMS", getBAMSURI() );
        int allenStructCount = 0;

        // get all BAMS regions, try to find match class
        for ( OntClass BAMSClass : BAMSData.getAllBrianRegions() ) {
            String BAMSClassLabel = BAMSClass.getLabel( null );

            // get the >200 structure expression mappings for this BAMS region (if any)
            // how many have more than one? use StructureCatalogLoader to test
            Set<String> allenCatalogNames = allenStructureLoader.getAllenMappedRegions( BAMSClassLabel );
            Set<String> allenCatalogExpressionNames = new HashSet<String>( catalogExpression.getEnergies()
                    .getColNames() );
            
            // get rid of the ones with no expression for all regions
            Set<String> removeZeroes = Util.findZeroColumns( catalogExpression.getEnergies() );
            allenCatalogExpressionNames.removeAll( Util.findZeroColumns( catalogExpression.getEnergies() ) );
            log.info( "removing all zero cols - " + removeZeroes );
            
            // log.info(allenCatalogExpressionNames);
            if ( allenCatalogNames != null ) {
                for ( String allenName : allenCatalogNames ) {
                    System.out.println( "allenStructureLoader HIT" );
                    Resource region = allenDataRDFModel.createResource( BAMSClass.getURI() );
                    allenStructCount++;
                    region.addProperty( hasStructureCatalogMapping, allenName );
                    if ( allenCatalogExpressionNames.contains( allenName ) ) {
                        System.out.println( "allenStructureDate HIT" );
                        region.addLiteral( Vocabulary.hasStructureCatalogExpression, true );
                    }
                }
            }

            if ( allenCatalogNames != null ) {
                for ( String allenName : allenCatalogNames ) {
                    System.out.println( "allenStructureLoader HIT" );
                    Resource region = allenDataRDFModel.createResource( BAMSClass.getURI() );
                    allenStructCount++;
                    region.addProperty( hasStructureCatalogMapping, allenName );
                }
            }

            // Tie into the Allen matrices for the >200 structure expression regions
            // if(allenStructureLoader.

            // get the top50+major17 allen mappings for this BAMS region (if any)
            Set<String> allenNames = BAMStoAllenMapping.get( BAMSClassLabel );
            if ( allenNames != null ) {
                for ( String allenName : allenNames ) {
                    Resource region = allenDataRDFModel.createResource( BAMSClass.getURI() );

                    // System.out.println( classLabel );
                    // get the Allen Genes for top50
                    Set<String> genes = allenLoader.getGenesForDongRegion( allenName );
                    if ( genes != null ) {
                        for ( String gene : genes ) {
                            region.addProperty( hasAllenEnrichedGene, gene );
                        }
                    }

                    // Tie into the Allen matrices for the 17 major regions
                    if ( allen17Names.contains( allenName ) ) {
                        System.out.println( "HIT" );
                        region.addProperty( hasAllen17Mapping, allenName );
                    }

                }
            }
        }

        allenDataRDFModel.write( new FileWriter( SetupParameters.config.getString( "abams.dataFolder" )
                + "RDFAllenToBAMSLinks.rdf" ) );
        // allenData.write( System.out );
        log.info( "Done" );
        log.info( "Allen to BAMS Links:" + allenStructCount );
    }
}
