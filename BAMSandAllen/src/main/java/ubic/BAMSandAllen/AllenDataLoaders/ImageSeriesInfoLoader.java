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
package ubic.BAMSandAllen.AllenDataLoaders;

import java.io.FileReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.SetupParameters;
import ubic.basecode.dataStructure.StringToStringSetMap;
import au.com.bytecode.opencsv.CSVReader;

public class ImageSeriesInfoLoader {
    private static Log log = LogFactory.getLog( ImageSeriesInfoLoader.class.getName() );

    String filename;
    boolean headerLine;
    char sep;
    int imageseriesid = 0;
    int geneNamePos = 2;

    int genesymbolPos = 3;
    int planePOS = 13;
    int entrezGeneID = 23;
    Map<String, String> imageIDtoGene, imageIDtoName;
    Map<String, Integer> imageIDtoEntrezID;
    StringToStringSetMap planeToImageID;

    public ImageSeriesInfoLoader() throws Exception {
        filename = SetupParameters.config.getString( "abams.allen.imageseriesCSV" );
        // filename = "/home/leon/temp/allen/2009Data/imageseries.csv";
        sep = ',';
        headerLine = true;
        imageIDtoGene = new HashMap<String, String>();
        imageIDtoName = new HashMap<String, String>();
        imageIDtoEntrezID = new HashMap<String, Integer>();
        planeToImageID = new StringToStringSetMap();
        init();
    }

    public void printEntrezIDStats() {
        log.info( "Image ID size:" + imageIDtoGene.size() );
        log.info( "Image IDs with entrezID:" + imageIDtoEntrezID.size() );
        log.info( "Unique EntrezIDs:" + new HashSet<Integer>( imageIDtoEntrezID.values() ).size() );
        log.info( "Unique GeneSymbols:" + new HashSet<String>( imageIDtoName.values() ).size() );

    }

    public void init() throws Exception {
        CSVReader reader = new CSVReader( new FileReader( filename ), sep );
        String[] line;
        // eat the header line
        int linecount = 0;
        if ( headerLine ) {
            linecount++;
            line = reader.readNext();
        }
        while ( ( line = reader.readNext() ) != null ) {
            linecount++;
            imageIDtoGene.put( line[imageseriesid], line[genesymbolPos] );
            imageIDtoName.put( line[imageseriesid], line[geneNamePos] );
            int entrezID = Integer.parseInt( line[entrezGeneID] );
            if ( entrezID != 0 ) imageIDtoEntrezID.put( line[imageseriesid], entrezID );
            String plane = line[planePOS].toLowerCase();
            planeToImageID.put( plane, line[imageseriesid] );
        }
        reader.close();
        log.info( "image id size:" + imageIDtoGene.size() );
        log.info( "gene names size:" + getAllGenes().size() );
        log.info( "lines (including header):" + linecount );
        log.info( "Plane counts, sagittal:" + planeToImageID.get( "sagittal" ).size() + " coronal:"
                + planeToImageID.get( "coronal" ).size() );

    }

    public Set<String> getAllGenes() {
        return new HashSet<String>( imageIDtoGene.values() );
    }

    public Set<String> getAllRowsWithCoronalGenes() {
        // ugly
        Set<String> result = new HashSet<String>();
        for ( String gene : imageIDtoGene.values() ) {
            if ( hasCoronalImageFromRowName( gene ) ) result.addAll( getRowsFromGene( gene ) );
        }
        return result;
    }

    /*
     * convert row labels like Sem3a[12345] to 12345
     */
    public static String getImageIDFromRowLabel( String rowName ) {
        Pattern p = Pattern.compile( ".*\\[(\\d+)\\]" );
        Matcher m = p.matcher( rowName );
        m.matches();
        return m.group( 1 );
    }

    public String getPlaneFromRowName( String rowName ) {
        String imageID = getImageIDFromRowLabel( rowName );
        // log.info( "imageID:" + imageID );

        return getPlaneFromImageID( imageID );
    }

    public static String getGeneNameFromRowName( String rowName ) {
        String geneName = rowName.replaceAll( "\\[.*\\]", "" );
        return geneName;
    }

    public boolean hasCoronalImageFromRowName( String rowName ) {
        String geneName = ImageSeriesInfoLoader.getGeneNameFromRowName( rowName );
        return isInCoronalSet( geneName );
    }

    private boolean isInCoronalSet( String geneName ) {
        Set<String> allImageSets = getRowsFromGene( geneName );
        for ( String otherImage : allImageSets ) {
            if ( getPlaneFromRowName( otherImage ).equals( "coronal" ) ) return true;
        }
        return false;
    }

    public String getPlaneFromImageID( String imageID ) {
        Set<String> result = planeToImageID.whereIs( imageID );
        // not found in the image loader, null gene name it seems
        if ( result.size() == 0 ) {
            return "unknown";
        }
        assert ( result.size() == 1 );
        // log.info( imageID );
        return result.iterator().next();
    }

    public Integer getNCBIIDFromRowName( String rowName ) {
        String imageID = getImageIDFromRowLabel( rowName );
        // log.info( "imageID:" + imageID );
        return imageIDtoEntrezID.get( imageID );
    }

    public String getSymbolFromNCBIID( String NCBIID ) {
        for ( String key : imageIDtoEntrezID.keySet() ) {
            if ( imageIDtoEntrezID.get( key ).toString().equals( NCBIID ) ) {
                return getGeneFromImageID( key );
            }
        }
        return null;
    }

    public String getNameFromImageID( String imageID ) {
        return imageIDtoName.get( imageID );
    }

    public Map<String, String> getMappings() {
        return imageIDtoGene;
    }

    /*
     * returns null if it can't map it
     */
    public String getGeneFromImageID( String imageid ) {
        return imageIDtoGene.get( imageid );
    }

    public Set<String> getRowsFromGene( String gene ) {
        Set<String> result = new HashSet<String>();
        for ( String key : imageIDtoGene.keySet() ) {
            String value = imageIDtoGene.get( key );
            if ( value.equals( gene ) ) result.add( value + "[" + key + "]" );
        }
        return result;
    }

    public static void main( String[] args ) throws Exception {
        ImageSeriesInfoLoader test = new ImageSeriesInfoLoader();
        log.info( getImageIDFromRowLabel( "rtest2[123456]" ) );
        log.info( test.getPlaneFromRowName( "Sema3g[69526474]" ) );
        test.printEntrezIDStats();
    }
}
