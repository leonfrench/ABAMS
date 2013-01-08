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

import static ubic.BAMSandAllen.Vocabulary.direct_part_of;
import static ubic.BAMSandAllen.Vocabulary.has_direct_part;

import java.io.FileInputStream;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CopyOnWriteArraySet;

import ubic.BAMSandAllen.JenaUtil;
import ubic.BAMSandAllen.LexiconSource;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Vocabulary;
import ubic.BAMSandAllen.AllenDataLoaders.RegionLoader;
import ubic.BAMSandAllen.AllenDataLoaders.StructureCatalogLoader;
import ubic.BAMSandAllen.ClassSelectors.BrainRegionClassSelector;
import ubic.BAMSandAllen.ClassSelectors.ClassSelector;

import com.hp.hpl.jena.ontology.AllValuesFromRestriction;
import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.ontology.OntDocumentManager;
import com.hp.hpl.jena.ontology.OntModel;
import com.hp.hpl.jena.ontology.OntModelSpec;
import com.hp.hpl.jena.ontology.Restriction;
import com.hp.hpl.jena.ontology.SomeValuesFromRestriction;
import com.hp.hpl.jena.rdf.model.Model;
import com.hp.hpl.jena.rdf.model.ModelFactory;
import com.hp.hpl.jena.rdf.model.NodeIterator;
import com.hp.hpl.jena.rdf.model.RDFNode;
import com.hp.hpl.jena.rdf.model.Resource;
import com.hp.hpl.jena.vocabulary.RDF;
import com.hp.hpl.jena.vocabulary.RDFS;

public class BAMSDataLoader extends RegionLoader implements LexiconSource {

    Set<OntClass> swansonRegions;
    OntModel BAMSModel;
    ClassSelector regionSelector;

    public OntModel getBAMSModel() {
        return BAMSModel;
    }

    public BAMSDataLoader() {
        try {
            OntDocumentManager mgr = OntDocumentManager.getInstance();
            mgr.setProcessImports( false );

            // set the mgr's properties now
            // now use it

            OntModelSpec s = new OntModelSpec( OntModelSpec.OWL_MEM );
            s.setDocumentManager( mgr );
            OntModel brainLinks = ModelFactory.createOntologyModel( s );

            brainLinks.read( new FileInputStream( SetupParameters.config.getString( "abams.brainLinksFile" ) ), "" );
            // brainLinks.read( Vocabulary.brainLinksURI );

            BAMSModel = ModelFactory.createOntologyModel( s );
            BAMSModel.read( new FileInputStream( SetupParameters.config.getString( "abams.BAMSDataInOWL" ) ),
                    Vocabulary.getBAMSURI() );

            // get all swanson regions, try to find match class
            regionSelector = new BrainRegionClassSelector();
            swansonRegions = JenaUtil.fitlerClasses( BAMSModel, regionSelector );
            // flatten the RDF
            makeFlat();
            // undo the process import change so it doesnt mess up the birnlex loader
            mgr.setProcessImports( true );
        } catch ( Exception e ) {
            e.printStackTrace();
            System.exit( 1 );
        }
    }

    public Set<OntClass> getAllBrianRegions() {
        return swansonRegions;
    }

    public Set<String> getRegionsForLexicon() {
        return getRegions();
    }

    public Set<String> getRegions() {
        Set<String> BAMSregions = new HashSet<String>();
        for ( OntClass region : getAllBrianRegions() ) {
            BAMSregions.add( region.getLabel( null ) );
        }
        return BAMSregions;
    }

    /*
     * populates a model for BAMS entries, copies information from BAMS ontology instead of linking and loading
     */
    public void addToModel( Model model ) {

        for ( OntClass region : getAllBrianRegions() ) {
            // any synonyms in BAMS? don't think there is any
            Resource r = model.createResource( region.getURI() );
            r.addLiteral( RDFS.label, region.getLabel( null ) );
            r.addProperty( RDF.type, Vocabulary.BAMSName );
            r.addProperty( Vocabulary.has_label_term, Vocabulary.makeNeurotermNode( region.getLabel( null ), model ) );
        }
    }

    /*
     * A bit ugly, converts to OntClass, then back to string for output (non-Javadoc)
     * 
     * @see ubic.BAMSandAllen.AllenDataLoaders.RegionLoader#getParent(java.lang.String)
     */
    public String getParent( String region ) {
        OntClass regionClass = convertStringToClassRegion( region );
        OntClass parent = getParent( regionClass );
        return convertClassRegionToString( parent );
    }

    public Set<String> getParents( String region ) {
        Set<String> result = getParentsRecursive( region );
        result.remove( region ); // dont count self
        return result;
    }

    public Set<String> getParentsRecursive( String region ) {
        if ( region == null ) return new HashSet<String>();
        String parent = getParent( region );
        Set<String> result;
        if ( parent != null ) {
            result = getParentsRecursive( parent );
        } else {
            result = new HashSet<String>();
        }
        result.add( region );
        return result;
    }

    public OntClass getParent( OntClass region ) {
        OntClass parent;
        RDFNode node = region.getPropertyValue( direct_part_of );
        if ( node != null ) {
            parent = ( OntClass ) ( node.as( OntClass.class ) );
        } else {
            parent = null;
        }
        return parent;
    }

    public Set<String> getChildren( String region, boolean indirect ) {
        Set<String> result = new HashSet<String>();
        Set<OntClass> children = getChildren( convertStringToClassRegion( region ), indirect );
        for ( OntClass regionClass : children ) {
            result.add( convertClassRegionToString( regionClass ) );
        }
        return result;
    }

    public Set<OntClass> getChildren( OntClass region ) {
        return getChildren( region, false );
    }

    public Set<String> getChildren( String region ) {
        return getChildren( region, false );
    }

    public Set<OntClass> getChildren( OntClass region, boolean inDirect ) {
        Set<OntClass> children = new HashSet<OntClass>();
        NodeIterator nodes = region.listPropertyValues( has_direct_part );
        // some of the has direct parts are collections, not regions, so check
        while ( nodes.hasNext() ) {
            RDFNode node = nodes.nextNode();
            // convert to ontclass
            OntClass BAMSClass = ( OntClass ) ( node.as( OntClass.class ) );
            if ( regionSelector.test( BAMSClass ) ) {
                children.add( BAMSClass );
            }

        }
        // a bit ugly here, recursion
        if ( inDirect ) {
            Set<OntClass> childrenChildren = new HashSet<OntClass>();
            childrenChildren.addAll( children );
            for ( OntClass child : children ) {
                childrenChildren.addAll( getChildren( child, true ) );
            }
            return childrenChildren;
        } else {

            return children;
        }
    }

    public OntClass convertStringToClassRegion( String region ) {
        for ( OntClass regionClass : swansonRegions ) {
            String label = regionClass.getLabel( null );
            if ( label.equals( region ) ) return regionClass;
        }
        return null;
    }

    public String convertClassRegionToString( OntClass region ) {
        if ( region == null ) return null;
        String result = region.getLabel( null );
        if ( result == null ) result = null;
        return result; // + "(" + ontClass.getLocalName() + ")";
    }

    public List<String> convertRegionstoNames( Set<OntClass> regions ) {
        Set<String> enrichedRegionNames = new CopyOnWriteArraySet<String>();
        for ( OntClass region : regions )
            enrichedRegionNames.add( convertClassRegionToString( region ) );
        List<String> list = new LinkedList<String>( enrichedRegionNames );
        Collections.sort( list );
        return list;
    }

    /*
     * turn the "some values from" and "all values from" subclasses into just RDF
     */
    public void makeFlat() {
        for ( OntClass ontClass : getAllBrianRegions() ) {
            // System.out.println( classToString( ontClass ) );

            for ( Object superClassobj : ontClass.listSuperClasses().toList() ) {

                OntClass superClass = ( OntClass ) superClassobj;
                // String superLabel = superClass.getLabel( null );

                if ( superClass.isRestriction() ) {
                    Restriction r = superClass.asRestriction();
                    // System.out.print( " -R>" + r.getOnProperty().getLocalName() + " -> " );
                    OntClass value;
                    if ( r.isAllValuesFromRestriction() ) {
                        AllValuesFromRestriction rAV = r.asAllValuesFromRestriction();
                        value = ( OntClass ) rAV.getAllValuesFrom();
                        // bad stuff here, just add it to the model
                        ontClass.addProperty( r.getOnProperty(), value );
                        // System.out.println( classToString( value ) );
                    } else if ( r.isSomeValuesFromRestriction() ) {
                        SomeValuesFromRestriction rSV = r.asSomeValuesFromRestriction();
                        value = ( OntClass ) rSV.getSomeValuesFrom();
                        // bad stuff here, just add it to the model
                        ontClass.addProperty( r.getOnProperty(), value );
                        // System.out.println( classToString( value ) );
                    } else {
                        // System.out.println( "?" );
                    }
                } else {
                    /*
                     * System.out.println( " ->" + classToString( superClass ) ); System.out.println( " ->" +
                     * superClass.getURI() ); if ( ontClass.hasSuperClass( Vocabulary.Mouse_Brain_Part ) )
                     * System.out.println( "TRUE SUPERCLASS" );
                     */
                }
            }
        }
    }

    public static void main( String args[] ) throws Exception {
        BAMSDataLoader bamsLoader = new BAMSDataLoader();
        Set<String> BAMSRegions = bamsLoader.getRegions();
        String test = BAMSRegions.iterator().next();
        System.out.println( test );
        System.out.println( bamsLoader.getParent( test ) );
        System.out.println( bamsLoader.getParents( test ) );

        int totalDepth = 0;
        for ( String region : BAMSRegions ) {
            totalDepth += bamsLoader.getParents( region ).size();
        }
        System.out.println( "Average depth:" + ( totalDepth / (double)BAMSRegions.size() ) );

        // StructureCatalogLoader dong = new StructureCatalogLoader();
        // Set<String> dongRegions = dong.getRegions();
        //
        // for ( String dongRegion : dongRegions ) {
        // for ( String bams : BAMSRegions ) {
        // if ( dongRegion.toLowerCase().trim().equals( bams.toLowerCase().trim() ) ) {
        // if ( dong.getBAMSMappedRegions( dongRegion ) == null ) {
        // System.out.println( bams );
        // }
        // }
        // }
        // }

    }
}
