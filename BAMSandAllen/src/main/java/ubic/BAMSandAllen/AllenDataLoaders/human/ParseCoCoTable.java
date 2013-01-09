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
 ******************************************************************************/package ubic.BAMSandAllen.AllenDataLoaders.human;

import java.net.URLEncoder;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

import ubic.BAMSandAllen.Vocabulary;
import ubic.basecode.util.FileTools;

import com.hp.hpl.jena.rdf.model.Model;
import com.hp.hpl.jena.rdf.model.Resource;
import com.hp.hpl.jena.vocabulary.RDF;
import com.hp.hpl.jena.vocabulary.RDFS;

public class ParseCoCoTable {
    Set<String> CoCoRegions;
    Map<String, String> CoCoSquares;

    /**
     * @param args
     */
    // change to union of regions
    public static void main( String[] args ) throws Exception {
        ParseCoCoTable table = new ParseCoCoTable();
        table.parseTable();
        System.out.println( "---------" );
        table.mapABATerms();
        System.out.println("STEMMING DISABLED");

    }

    public void parseTable() throws Exception {
        // TODO Auto-generated method stub
        CoCoSquares = new HashMap<String, String>();
        CoCoRegions = new HashSet<String>();
        List<String> lines = FileTools.getLines( "C:\\Users\\leon\\Desktop\\CoCo\\tableS1Edit.txt" );
        for ( String line : lines ) {
            StringTokenizer tokes = new StringTokenizer( line, " ;" );
            LinkedList<String> tokens = new LinkedList<String>();
            while ( tokes.hasMoreElements() ) {
                tokens.addLast( tokes.nextToken().trim() );
            }
            int end = -1;
            for ( int i = 0; i < tokens.size() - 2; i++ ) {
                try {
                    if ( tokens.get( i + 3 ).contains( "[" ) ) break;
                    Integer.parseInt( tokens.get( i ) );
                    Integer.parseInt( tokens.get( i + 1 ) );
                    Integer.parseInt( tokens.get( i + 2 ) );
                    Integer.parseInt( tokens.get( i + 3 ) );
                    end = i;
                    // break;
                } catch ( Exception e ) {

                }
            }
            if ( end == -1 ) {
                System.out.println( "Failed on:" + line );
                System.out.println( tokens );
            }
            String outputLine = "\"";
            String regionName = "";
            for ( int i = 2; i < end; i++ ) {
                String token = tokens.get( i );

                regionName += token + " ";
            }
            regionName = regionName.trim();
            outputLine += regionName.trim();
            outputLine += "\"";
            // System.out.println( outputLine );
            int level = Integer.parseInt( tokens.get( 0 ).trim() );
            String acro = tokens.get( 1 ).trim();
            int ring = Integer.parseInt( tokens.get( end ).trim() );
            int degree = Integer.parseInt( tokens.get( end + 1 ).trim() );
            int timesStudied = Integer.parseInt( tokens.get( end + 2 ).trim() );
            int numberOfReportedConnections = Integer.parseInt( tokens.get( end + 3 ).trim() );

            CoCoRegions.add( regionName );

            String squared = "";
            for ( int i = end + 4; i < tokens.size(); i++ ) {
                String token = tokens.get( i );

                squared += token + " ";
            }
            CoCoSquares.put( regionName, squared );

            // System.out.println(squared);

        }

    }

    public void mapABATerms() throws Exception {
        List<String> ABARegions = FileTools.getLines( "C:\\Users\\leon\\Desktop\\Human array\\unionRegions.txt" );
        int hitCount = 0;
        Set<String> usedCoCoRegions = new HashSet<String>();
        for ( String ABARegion : ABARegions ) {
            String ABAOriginal = ABARegion;
            ABARegion = ABARegion.trim();
            ABARegion = ABARegion.replace( "Left", "" );
            ABARegion = ABARegion.replace( "Right", "" );

            Set<String> ABABag = makeBagOfWords( ABARegion );

            boolean hit = false;
            for ( String CoCoRegion : CoCoRegions ) {
                Set<String> CoCoBag = makeBagOfWords( CoCoRegion );
                if ( CoCoBag.equals( ABABag ) ) {
                    hit = true;
                    System.out.println( ABAOriginal + "\t" + CoCoRegion + "\t" + CoCoSquares.get( CoCoRegion ) );
                    usedCoCoRegions.add( CoCoRegion );
                }
            }
            if ( hit == false ) {
                System.out.println( ABAOriginal );
            } else {
                hitCount++;
            }

        }
        Set<String> unusedCoCo = new HashSet<String>( CoCoRegions );
        unusedCoCo.removeAll( usedCoCoRegions );
        for ( String CoCoRegion : unusedCoCo ) {
            System.out.println( "\t" + CoCoRegion + "\t" + CoCoSquares.get( CoCoRegion ) );
        }
        System.out.println( "Hitcount:" + hitCount );
    }

    public void addABATerms( Model model ) throws Exception {
        List<String> ABARegions = FileTools.getLines( "C:\\Users\\leon\\Desktop\\Human array\\unionRegions.txt" );

        Set<String> ABASet = new HashSet<String>();
        for ( String ABARegion : ABARegions ) {
            // ABARegion = ABARegion.replace( "Left", "" );
            // ABARegion = ABARegion.replace( "Right", "" );
            ABARegion = ABARegion.trim();
            ABASet.add( ABARegion );
        }

        String URIBase = "http://allenbrainatlas.com/humanRegions#";
        for ( String ABARegion : ABASet ) {
            String uri = URIBase + URLEncoder.encode( ABARegion, "UTF-8" );
            String label = ABARegion;
            Resource mainConcept = model.createResource( uri );
            mainConcept.addProperty( RDF.type, Vocabulary.ABAName );
            Resource r = Vocabulary.makeNeurotermNode( label, model );
            mainConcept.addLiteral( RDFS.label, label );
            mainConcept.addProperty( Vocabulary.has_label_term, r );
            System.out.println( uri + " = " + ABARegion );
        }
        System.out.println( "Regions:" + ABASet.size() );

    }

    // ugly
    // static Stemmer stemmer;
    // static {
    // stemmer = new LovinsStemmer();
    // }
    public static String delims = "~`!@#$%^&*()_-+={[}]|\\:;\"',<.>?/ \t\n\r";

    public Set<String> makeBagOfWords( String phrase ) {

        phrase = phrase.toLowerCase().trim();
        // phrase = stemmer.stemString( phrase );
        StringTokenizer tokens = new StringTokenizer( phrase, delims, false );
        Set<String> bagOfWords = new HashSet<String>();
        while ( tokens.hasMoreTokens() ) {
            bagOfWords.add( tokens.nextToken() );
        }
        bagOfWords.remove( "th" );
        bagOfWords.remove( "of" );
        return bagOfWords;
    }

}
