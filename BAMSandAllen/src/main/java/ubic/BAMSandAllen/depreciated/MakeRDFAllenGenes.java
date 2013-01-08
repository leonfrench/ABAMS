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
package ubic.BAMSandAllen.depreciated;

import static ubic.BAMSandAllen.Vocabulary.getBAMSURI;
import static ubic.BAMSandAllen.Vocabulary.getBrainLinksURI;
import static ubic.BAMSandAllen.Vocabulary.hasAllenEnrichedGene;

import java.io.FileWriter;
import java.util.Map;
import java.util.Set;

import ubic.BAMSandAllen.FullMappingLoader;
import ubic.BAMSandAllen.AllenDataLoaders.AllenTop50DataLoader;
import ubic.BAMSandAllen.BAMSDataLoaders.BAMSDataLoader;

import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.rdf.model.Model;
import com.hp.hpl.jena.rdf.model.ModelFactory;
import com.hp.hpl.jena.rdf.model.Resource;

public class MakeRDFAllenGenes {

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // load the mappings - CSV

        Map<String, Set<String>> BAMStoAllenMapping = new FullMappingLoader().getAllentoBamsMapping();

        BAMSDataLoader BAMSData = new BAMSDataLoader();

        Model allenData = ModelFactory.createDefaultModel();

        allenData.setNsPrefix( "brianLink", getBrainLinksURI() );
        allenData.setNsPrefix( "BAMS", getBAMSURI() );

        // get all BAMS regions, try to find match class
        AllenTop50DataLoader allenLoader = new AllenTop50DataLoader();

        for ( OntClass ontClass : BAMSData.getAllBrianRegions() ) {
            String classLabel = ontClass.getLabel( null );
            Set<String> allenNames = BAMStoAllenMapping.get( classLabel );
            if ( allenNames != null ) {
                for ( String allenName : allenNames ) {
                    // System.out.println( classLabel );
                    // get the Allen Genes
                    Set<String> genes = allenLoader.getGenesForDongRegion( allenName );
                    Resource region = allenData.createResource( ontClass.getURI() );

                    for ( String gene : genes ) {
                        region.addProperty( hasAllenEnrichedGene, gene );
                    }
                }
            }
        }

        allenData.write( new FileWriter( "data\\RDFAllenGenes.rdf" ) );
        allenData.write( System.out );
    }
}
