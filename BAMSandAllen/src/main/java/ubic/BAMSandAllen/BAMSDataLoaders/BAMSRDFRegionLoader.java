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
package ubic.BAMSandAllen.BAMSDataLoaders;

import java.io.File;
import java.io.FileInputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.JenaUtil;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.BAMSandAllen.AllenDataLoaders.RegionLoader;
import ubic.basecode.dataStructure.CountingMap;
import ubic.basecode.dataStructure.StringToStringSetMap;

import com.hp.hpl.jena.rdf.model.Model;
import com.hp.hpl.jena.rdf.model.ModelFactory;
import com.hp.hpl.jena.rdf.model.Property;
import com.hp.hpl.jena.rdf.model.Resource;
import com.hp.hpl.jena.vocabulary.OWL;
import com.hp.hpl.jena.vocabulary.RDF;

public class BAMSRDFRegionLoader {
    protected static Log log = LogFactory.getLog( BAMSRDFRegionLoader.class );
    Model model;

    Property class1;
    Property class2;
    private StringToStringSetMap parents;

    public BAMSRDFRegionLoader() throws Exception {
        parents = new StringToStringSetMap();
        String filename = SetupParameters.config.getString( "whitetext.lexicon.output" ) + "bams-serialization.xml";
        model = ModelFactory.createDefaultModel();
        File folder = new File( SetupParameters.config.getString( "abams.NIF.BAMS.sparql" ) );
        model.read( new FileInputStream( filename ), null );
        log.info( "Model size:" + model.size() );
        Property name = model.createProperty( "http://brancusi1.usc.edu/RDF/name" );
        Property species = model.createProperty( "http://brancusi1.usc.edu/RDF/species" );
        class1 = model.createProperty( "http://brancusi1.usc.edu/RDF/class1" );
        class2 = model.createProperty( "http://brancusi1.usc.edu/RDF/class2" );
        Resource partType = model.createResource( "http://brancusi1.usc.edu/RDF/brainPart" );

        Resource hierType = model.createResource( "http://brancusi1.usc.edu/RDF/hierarchy" );
        Set<Resource> parts = JenaUtil.getSubjects( model.listStatements( null, RDF.type, partType ) );
        Set<Resource> hiers = JenaUtil.getSubjects( model.listStatements( null, RDF.type, hierType ) );
        Set<String> regions = new HashSet<String>();
        for ( Resource part : parts ) {
            String label = JenaUtil.getStringLiteral( part, name );
            // if ( label != null )
            String speciesString = JenaUtil.getStringLiteral( part, species );
            // log.info( label );
            if ( speciesString.equals( "Rat" ) ) {
                regions.add( label );
            }
        }

        // <bams:class1 rdf:resource="http://brancusi1.usc.edu/brain_parts/nucleus-of-the-solitary-tract-6/"/>
        // <bams:class2
        // rdf:resource="http://brancusi1.usc.edu/brain_parts/nucleus-of-the-solitary-tract-commissural-part-5/"/>
        // <rdf:type rdf:resource="http://brancusi1.usc.edu/RDF/hierarchy"/>
        // <owl:OntologyProperty rdf:resource="http://brancusi1.usc.edu/RDF/hasPart"/>
        // check organism!! and existence of region??
        for ( Resource heirPair : hiers ) {
            Resource class1Res = heirPair.getProperty( class1 ).getObject().as( Resource.class );
            Resource class2Res = heirPair.getProperty( class2 ).getObject().as( Resource.class );
            String class1String = JenaUtil.getStringLiteral( class1Res, name );
            String class2String = JenaUtil.getStringLiteral( class2Res, name );
            String class1Species = JenaUtil.getStringLiteral( class1Res, species );
            String class2Species = JenaUtil.getStringLiteral( class2Res, species );

            if ( class2Species.equals( "Rat" ) && class1Species.equals( "Rat" ) ) {
                parents.put( class2String, class1String );

                String relationship = heirPair.getProperty( OWL.OntologyProperty.as( Property.class ) ).getObject()
                        .toString();
                if ( !relationship.endsWith( "hasPart" ) ) {
                    log.info( relationship );
                }
            }
        }
    }

    public void fontiersTesting() {
        Resource topRelationType = model.createResource( "http://brancusi1.usc.edu/RDF/topologicalRelation" );
        Property nomenclature = model.createProperty( "http://brancusi1.usc.edu/RDF/nomenclature" );
        Property name = model.createProperty( "http://brancusi1.usc.edu/RDF/name" );
        Resource partType = model.createResource( "http://brancusi1.usc.edu/RDF/brainPart" );

        Set<Resource> topRelations = JenaUtil.getSubjects( model.listStatements( null, RDF.type, topRelationType ) );
        CountingMap<String> counts = new CountingMap<String>();

        // get 98 regions
        // get 2004 regions
        Set<Resource> parts = JenaUtil.getSubjects( model.listStatements( null, RDF.type, partType ) );
        Set<Resource> regions98 = new HashSet<Resource>();
        Set<Resource> regions04 = new HashSet<Resource>();
        for ( Resource part : parts ) {
            Resource nomenRes = part.getProperty( nomenclature ).getObject().as( Resource.class );
            String nomenName = JenaUtil.getStringLiteral( nomenRes, name );
            if ( nomenName.equals( "Swanson-1998" ) ) {
                regions98.add( part );
            }
            if ( nomenName.equals( "Swanson-2004" ) ) {
                regions04.add( part );
            }
        }

        Set<Resource> seen98 = new HashSet<Resource>();
        Set<Resource> seen04 = new HashSet<Resource>();
        Set<String> strings98 = new HashSet<String>();
        Set<String> strings04 = new HashSet<String>();

        for ( Resource topRelation : topRelations ) {
            Resource class1Res = topRelation.getProperty( class1 ).getObject().as( Resource.class );
            Resource class2Res = topRelation.getProperty( class2 ).getObject().as( Resource.class );
            Resource relationship = topRelation.getProperty( OWL.OntologyProperty.as( Property.class ) ).getObject()
                    .as( Resource.class );
            Resource nomen1 = class1Res.getProperty( nomenclature ).getObject().as( Resource.class );
            Resource nomen2 = class2Res.getProperty( nomenclature ).getObject().as( Resource.class );
            log.info( nomen1.toString() );
            log.info( nomen2.toString() );
            log.info( JenaUtil.getObjects( nomen1.listProperties() ).toString() );
            log.info( JenaUtil.getSubjects( nomen1.listProperties() ).toString() );
            String nomenName1 = JenaUtil.getStringLiteral( nomen1, name );
            String nomenName2 = JenaUtil.getStringLiteral( nomen2, name );

            String class1String = JenaUtil.getStringLiteral( class1Res, name );
            String class2String = JenaUtil.getStringLiteral( class2Res, name );

            // filter for swanson 1998 to 2004
            Set<String> nomens = new HashSet<String>();
            nomens.add( nomenName1 );
            nomens.add( nomenName2 );

            if ( nomens.contains( "Swanson-1998" ) && nomens.contains( "Swanson-2004" ) ) {
                // quantify <owl:OntologyProperty
                counts.increment( relationship.toString() );

                if ( nomenName1.equals( "Swanson-1998" ) ) {
                    seen98.add( class1Res );
                    strings98.add( class1String.toLowerCase() );
                    seen04.add( class2Res );
                    strings04.add( class2String.toLowerCase() );
                } else {
                    seen98.add( class2Res );
                    strings98.add( class2String.toLowerCase() );
                    seen04.add( class1Res );
                    strings04.add( class1String.toLowerCase() );

                }

            }
        }
        log.info( counts.toString() );

        log.info( "Seen 1998 mappings:" + seen98.size() );
        log.info( "Seen 2004 mappings:" + seen04.size() );
        log.info( "Total 1998:" + regions98.size() );
        log.info( "Total 2004:" + regions04.size() );

        log.info( "Total strings 1998:" + strings98.size() );
        log.info( "Total strings 2004:" + strings04.size() );

        log.info( "string intersect:" + Util.intersectSize( strings98, strings04 ) );

        // log.info( "Unseen 1998:"+ Util.subtract( regions98, seen98 ) );

    }

    public Set<String> getDirectChildren( String region ) {
        Set<String> result = new HashSet<String>();
        for ( String possibleChild : parents.keySet() ) {
            if ( parents.get( possibleChild ).contains( region ) ) {
                result.add( possibleChild );
            }
        }
        return result;
    }

    public Set<String> getParent( String region ) {
        return parents.get( region );
    }

    public Set<String> getRegions() {
        Set<String> result = new HashSet<String>();
        result.addAll( parents.getSeenValues() );
        result.addAll( parents.keySet() );
        return result;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        BAMSRDFRegionLoader loader = new BAMSRDFRegionLoader();
        log.info( loader.getParent( "Substantia nigra" ) );
        log.info( loader.getParent( "Behavioral state system" ) );
        log.info( loader.getDirectChildren( "Behavioral state system" ) );
        log.info( loader.getDirectChildren( "Substantia nigra" ) );
        log.info( loader.getDirectChildren( "Substantia nigra compact part" ) );
        log.info( loader.getParent( "Cerebrospinal trunk" ) );
        log.info( loader.getParent( "Central nervous system gray matter" ) );
        log.info( loader.getParent( "Central Nervous System" ) == null );
        loader.fontiersTesting();

    }

}
