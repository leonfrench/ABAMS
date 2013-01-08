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
package ubic.BAMSandAllen.MappingHelpers;

import java.io.FileInputStream;
import java.util.HashSet;
import java.util.Set;


import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.ontology.OntDocumentManager;
import com.hp.hpl.jena.ontology.OntModel;
import com.hp.hpl.jena.ontology.OntModelSpec;
import com.hp.hpl.jena.rdf.model.ModelFactory;
import com.hp.hpl.jena.util.iterator.ExtendedIterator;

public class JenaTry {
    public static void main( String args[] ) throws Exception {
        OntDocumentManager mgr = new OntDocumentManager();
        // set the mgr's properties now
        // now use it
        OntModelSpec s = new OntModelSpec( OntModelSpec.OWL_MEM );
        s.setDocumentManager( mgr );
        OntModel brainLinks = ModelFactory.createOntologyModel( s );
        brainLinks.read(
                new FileInputStream( "C:\\Documents and Settings\\lfrench\\Desktop\\ontoloties\\brainLinks.owl" ), "" );

        OntModel bams = ModelFactory.createOntologyModel( s );
        bams.read( new FileInputStream(
                "C:\\Documents and Settings\\lfrench\\Desktop\\ontoloties\\bams-from-swanson-98-4-23-07.owl" ), "" );

        // BamsAllenMapper fullMapper = new BamsAllenMapper();

        Set<String> swansonRegions = BamsAllenMapper.getAllSwansonNames();
        Set<String> seen = new HashSet<String>();

        // get all swanson regions, try to find match class
        ExtendedIterator ei = bams.listClasses();
        int count = 0;
        for ( Object o = ei.next(); ei.hasNext(); o = ei.next() ) {
            OntClass ontClass = ( OntClass ) o;
            String classLabel = ontClass.getLabel( null );
            if ( classLabel != null ) {
                // System.out.println(classLabel);
                boolean found = false;
                for ( String swansonRegion : swansonRegions ) {
                    if ( swansonRegion.replaceAll( ",", "" ).equalsIgnoreCase( classLabel.replaceAll( ",", "" ) ) ) {
                        count++;
                        seen.add( swansonRegion );
                        //System.out.println( classLabel );
                        System.out.println( swansonRegion +"->"+ classLabel );
                        found = true;
                    }
                }
            }
        }
        System.out.println( count + " of " + swansonRegions.size()
                + " Allen enriched regions were found in the BAMS dataset" );
        System.out.println( "not found:" );
        swansonRegions.removeAll( seen );
        for ( String ss : swansonRegions ) {
            System.out.println( ss );
        }

    }
}
