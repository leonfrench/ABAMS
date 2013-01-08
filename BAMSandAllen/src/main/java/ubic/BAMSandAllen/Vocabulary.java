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

import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;

import com.hp.hpl.jena.ontology.OntClass;
import com.hp.hpl.jena.rdf.model.Model;
import com.hp.hpl.jena.rdf.model.ModelFactory;
import com.hp.hpl.jena.rdf.model.Property;
import com.hp.hpl.jena.rdf.model.Resource;
import com.hp.hpl.jena.vocabulary.RDF;
import com.hp.hpl.jena.vocabulary.RDFS;

public class Vocabulary {
    private static Model m;
    static {
        m = ModelFactory.createDefaultModel();
    }

    protected static final String BAMSURI = "http://hissa.nist.gov/jb/biordf-demo/bams-from-swanson-98-4-23-07.owl#";

    protected static final String ABAURI = "http://community.brain-map.org/confluence/download/attachments/525267/ABA.owl?version=1#";

    protected static final String brainLinksURI = "http://ubic.ubc.ca/ontology/brainLinks.owl#";

    protected static final String lexiconLinksURI = "http://www.chibi.ubc.ca/Gemma/ws/xml/neuroanatomyLinks.owl#";

    protected static final String NNURI = "http://www.purl.org/neuronames#";

    protected static final String BredeURI = "http://neuro.imm.dtu.dk/services/brededatabase/worois.xml#";

    protected static final String lexiconSpaceURI = "http://www.purl.org/lexiconspace#";

    protected static final String mentionSpaceURI = "http://www.purl.org/mentionspace#";

    protected static final String pubmedURIPrefix = "http://bio2rdf.org/pubmed:";

    public static String getpubmedURIPrefix() {
        return pubmedURIPrefix;
    }

    public static String getBAMSURI() {
        return BAMSURI;
    }

    public static String getMentionSpaceURI() {
        return mentionSpaceURI;
    }

    public static String getABAURI() {
        return ABAURI;
    }

    public static String getLexiconURI() {
        return lexiconLinksURI;
    }

    public static String getNNURI() {
        return NNURI;
    }

    public static String getBredeURI() {
        return BredeURI;
    }

    public static String getBrainLinksURI() {
        return brainLinksURI;
    }

    public static String getLexiconSpaceURI() {
        return lexiconSpaceURI;
    }

    public static boolean isEnriched( OntClass ontClass ) {
        return ontClass.hasProperty( hasAllenEnrichedGene );
    }

    public static boolean isTop17( OntClass ontClass ) {
        return ontClass.hasProperty( hasAllen17Mapping );
    }

    public static Resource makeNeurotermNode( String s, Model model ) {
        boolean toLower = true;
        return makeNeurotermNode( s, model, lexiconSpaceURI, Vocabulary.neuroterm, true, toLower );
    }

    public static Resource getNeurotermNode( String s, Model model ) {
        boolean toLower = true;
        return makeNeurotermNode( s, model, lexiconSpaceURI, Vocabulary.neuroterm, false, toLower );
    }

    public static Resource makeMentionNode( String s, Model model ) {
        boolean toLower = false;
        boolean newNode = true;
        return makeNeurotermNode( s, model, mentionSpaceURI, Vocabulary.neuromention, newNode, toLower );
    }

    public static Resource getMentionNode( String s, Model model ) {
        boolean toLower = false;
        boolean newNode = false;
        return makeNeurotermNode( s, model, mentionSpaceURI, Vocabulary.neuromention, newNode, toLower );
    }

    private static Resource makeNeurotermNode( String s, Model model, String URIBase, Resource type, boolean newNode,
            boolean toLower ) {
        Resource r = null;
        r = model.createResource( getMentionURI( s, URIBase, toLower ) );
        if ( newNode ) {
            r.addProperty( RDFS.label, s );
            // set as neuroterm class
            r.addProperty( RDF.type, type );
        }
        return r;
    }

    public static String getMentionURI( String s ) {
        boolean toLower = true;
        return getMentionURI( s, mentionSpaceURI, toLower );
    }

    public static String getMentionURI( String s, String URIBase, boolean toLower ) {
        String uri = null;
        s = s.trim();
        // convert to lowercase?
        if ( toLower ) s = s.toLowerCase();
        try {
            uri = URIBase + URLEncoder.encode( s, "UTF-8" );
        } catch ( UnsupportedEncodingException e ) {
            e.printStackTrace();
        }
        return uri;
    }

    // I don't know why it is not rat brain part
    public static final Resource Mouse_Brain_Part = m.createResource( BAMSURI + "Mouse_Brain_Part" );

    // may need to update
    public static final Resource Bed_Nuclei_Stria_Part = m.createResource( BAMSURI + "p104" );


    public static final Property direct_part_of = m.createProperty( BAMSURI + "direct_part_of" );
    public static final Property has_direct_part = m.createProperty( BAMSURI + "has_direct_part" );

    public static final Property hasAllenEnrichedGene = m.createProperty( brainLinksURI + "hasAllenEnrichedGene" );
    public static final Property hasAllen17Mapping = m.createProperty( brainLinksURI + "hasAllen17Mapping" );
    public static final Property hasStructureCatalogMapping = m.createProperty( brainLinksURI
            + "hasStructureCatalogMapping" );
    public static final Property hasStructureCatalogExpression = m.createProperty( brainLinksURI
            + "hasStructureCatalogExpression" );
    public static final Property has_direct_target = m.createProperty( BAMSURI + "has_direct_target" );
    public static final Property has_direct_source = m.createProperty( BAMSURI + "has_direct_source" );

    public static final Property mentions_species = m.createProperty( lexiconLinksURI + "mentions_species" );
    public static final Property publication_date = m.createProperty( lexiconLinksURI + "publication_date" );
    public static final Property has_NN_latin_term = m.createProperty( lexiconLinksURI + "has_NN_latin_term" );
    public static final Property has_label_term = m.createProperty( lexiconLinksURI + "has_label_term" );
    public static final Property has_synonym_term = m.createProperty( lexiconLinksURI + "has_synonym_term" );

    public static final Property has_Dong_term = m.createProperty( lexiconLinksURI + "has_Dong_term" );
    public static final Property has_Paxinos_Franklin_term = m.createProperty( lexiconLinksURI
            + "has_Paxinos_Franklin_term" );
    public static final Property has_Hof_term = m.createProperty( lexiconLinksURI + "has_Hof_term" );
    public static final Property has_Swanson_term = m.createProperty( lexiconLinksURI + "has_Swanson_term" );
    public static final Property has_NN_link = m.createProperty( lexiconLinksURI + "has_NN_link" );

    public static final Property in_PMID = m.createProperty( lexiconLinksURI + "in_PMID" );

    public static final Property has_manual_link = m.createProperty( lexiconLinksURI + "has_manual_link" );

    // from evaluations
    public static final Property evaluation_accept = m.createProperty( lexiconLinksURI + "evaluation_accept" );
    public static final Property evaluation_reject = m.createProperty( lexiconLinksURI + "evaluation_reject" );
    public static final Property evaluation_result = m.createProperty( lexiconLinksURI + "evaluation_result" );
    public static final Property evaluation_specific_to_general = m.createProperty( lexiconLinksURI
            + "evaluation_specific_to_general" );

    public static final Property number_of_occurances = m.createProperty( lexiconLinksURI + "number_of_occurances" );
    public static final Property number_of_abstracts = m.createProperty( lexiconLinksURI + "number_of_abstracts" );
    public static final Property annotation_set = m.createProperty( lexiconLinksURI + "annotation_set" );
    public static final Property annotated_connection_partner = m.createProperty( lexiconLinksURI
            + "annotated_connection_partner" );

    // match types
    public static final Property match = m.createProperty( lexiconLinksURI + "match" );

    public static final Property string_match_ignorecase = m.createProperty( lexiconLinksURI
            + "string_match_ignorecase" );
    public static final Property stem_match = m.createProperty( lexiconLinksURI + "stem_match" );
    public static final Property three_match = m.createProperty( lexiconLinksURI + "three_match" );
    public static final Property simple_mapping_match = m.createProperty( lexiconLinksURI + "simple_mapping_match" );
    public static final Property word_bag_match_ignorecase = m.createProperty( lexiconLinksURI
            + "word_bag_match_ignorecase" );
    public static final Property stem_bag_match_ignorecase = m.createProperty( lexiconLinksURI
            + "stem_bag_match_ignorecase" );

    // set as class
    public static final Resource neuroOntologyEntry = m.createResource( lexiconLinksURI + "neuroOntologyEntry" );
    public static final Resource neuroname = m.createResource( lexiconLinksURI + "neuroname" );
    public static final Resource bredeName = m.createResource( lexiconLinksURI + "bredeName" );
    public static final Resource ABAName = m.createResource( lexiconLinksURI + "ABAName" );
    // avian brain connectivity databsae
    public static final Resource ABCDName = m.createResource( lexiconLinksURI + "ABCDName" );
    public static final Resource BAMSName = m.createResource( lexiconLinksURI + "BAMSName" );
    public static final Resource neuromention = m.createResource( lexiconLinksURI + "neuromention" );
    public static final Resource BIRNLexname = m.createResource( lexiconLinksURI + "BIRNLexname" );
    public static final Resource neuroterm = m.createResource( lexiconLinksURI + "neuroterm" );
    
   

}
