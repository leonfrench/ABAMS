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

import java.io.Serializable;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.apache.poi.hssf.usermodel.HSSFSheet;

import ubic.BAMSandAllen.LexiconSource;
import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.BAMSDataLoaders.BAMSDataLoader;
import ubic.basecode.dataStructure.StringToStringSetMap;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import ubic.basecode.io.excel.ExcelUtil;

public class StructureCatalogLoader extends RegionLoader implements LexiconSource, Serializable {
    protected static Log log = LogFactory.getLog( StructureCatalogLoader.class );

    public static final int labelPosition = 0;
    public static final int namePosition = 1;
    public static final int parentLabelPosition = 2;
    public static final int mappedBAMSPosition = 3;
    public static final int motorNucleiPosition = 5;

    private Map<String, String> acronyms;
    private Map<String, String> parents;
    private StringToStringSetMap BAMSMap;
    private Set<String> motorNuclei;

    public StructureCatalogLoader() throws Exception {
        acronyms = new HashMap<String, String>();
        parents = new HashMap<String, String>();
        motorNuclei = new HashSet<String>();
        // a map<String,Set<String>>

        BAMSMap = new StringToStringSetMap();

        HSSFSheet sheet = ExcelUtil.getSheetFromFile( SetupParameters.config.getString( "abams.allenCDMapping" ),
                "Sheet1" );

        // eat the header line
        int row = 0;

        while ( true ) {
            row++; // skips the header line

            String label = ExcelUtil.getValue( sheet, row, labelPosition );
            String fullName = ExcelUtil.getValue( sheet, row, namePosition );
            String parentLabel = ExcelUtil.getValue( sheet, row, parentLabelPosition );
            String mappedBAMSRegion = ExcelUtil.getValue( sheet, row, mappedBAMSPosition );
            String motorRegion = ExcelUtil.getValue( sheet, row, motorNucleiPosition );

            if ( label == null ) break;

            if ( fullName != null ) fullName = fullName.trim();
            if ( mappedBAMSRegion != null ) mappedBAMSRegion = mappedBAMSRegion.trim();

            acronyms.put( label, fullName );
            if ( parentLabel != null ) {
                parents.put( fullName, parentLabel );
            }

            if ( mappedBAMSRegion != null && !mappedBAMSRegion.equals( "" ) ) {
                BAMSMap.put( fullName, mappedBAMSRegion.trim() );
            }

            if ( motorRegion.equals( "1.0" ) ) {
                motorNuclei.add( fullName );
            }
        }

        // convert parent acros to full forms
        for ( String key : parents.keySet() ) {
            String expandedParent = acronyms.get( parents.get( key ) );
            parents.put( key, expandedParent );
        }

    }

    /**
     * Given a BAMS region return the mapped allen regions
     */
    public Set<String> getAllenMappedRegions( String BAMSregion ) {
        Set<String> result = new HashSet<String>();
        for ( String allenRegion : BAMSMap.keySet() ) {
            if ( BAMSMap.get( allenRegion ).contains( BAMSregion ) ) result.add( allenRegion );
        }
        return result;
    }

    public void doStats() {
        System.out.println( "regions:" + getRegions().size() );
        System.out.println( "leafs:" + getLeafAmount() );
        System.out.println( "BAMS regions mapped to a Allen region:" + getBAMSMappedRegions().size() );
        System.out.println( "Allen regions that are mapped to a BAMS region:" + getAllenMappedRegions().size() );

        BAMSDataLoader bams = new BAMSDataLoader();

        Set<String> leftBAMS = bams.getRegions();
        System.out.println( "All BAMS regions:" + leftBAMS.size() );

        leftBAMS.removeAll( getBAMSMappedRegions() );
        System.out.println( "Remaining BAMS regions that are not mapped:" + leftBAMS.size() );

        int allenToMany = 0;
        int BAMSToMany = 0;
        for ( String allenRegion : BAMSMap.keySet() ) {
            if ( getBAMSMappedRegions( allenRegion ).size() > 1 ) {
                System.out.println( allenRegion );
                allenToMany++;
            }
        }
        System.out.println( "Allen regions pointing to more than one BAMS region:" + allenToMany );

        for ( String BAMSRegion : BAMSMap.getSeenValues() ) {
            if ( getAllenMappedRegions( BAMSRegion ).size() > 1 ) BAMSToMany++;
        }
        System.out.println( "BAMS regions pointing to more than one allen region:" + BAMSToMany );
    }

 

    public Set<String> getBAMSMappedRegions( String name ) {
        Set<String> result = BAMSMap.get( name );
        // make a new set incase it is modified
        if ( result != null ) result = new HashSet<String>( BAMSMap.get( name ) );
        return result;
    }

    /*
     * Case senstive retrivial
     */
    public String getFullName( String acronym ) {
        String fullName = acronyms.get( acronym );
        // special case for dupes
        if ( fullName.equals( "NO NAME" ) ) return acronym;
        if ( fullName.equals( "???" ) ) return acronym;
        return fullName;
    }

    /*
     * Case senstive retrivial
     */
    public String getAcro( String name ) {
        for ( String acro : acronyms.keySet() ) {
            if ( acronyms.get( acro ).equals( name ) ) return acro;
        }
        return null;
    }

    /*
     * Case senstive retrivial
     */
    public String getParent( String name ) {
        return parents.get( name );
    }

    public Set<String> getParents( String name ) {
        Set<String> result = new HashSet<String>();
        String parent = name;
        while ( null != ( parent = getParent( parent ) ) ) {
            result.add( parent );
        }
        return result;
    }

    public Set<String> getChildren( String name ) {
        Set<String> result = new HashSet<String>();
        // a problem occurs with Ammon's horn, as in terms of regions with data it is a leaf, but in the catalog loader
        // it is not, so return no children
        if ( name.equals( "Ammon\'s Horn" ) ) return result;
        for ( String key : parents.keySet() ) {
            String parent = parents.get( key );
            if ( parent != null && parent.equals( name ) ) {
                result.add( key );
            }
        }
        return result;
    }

    public Set<String> getRegionsForLexicon() {
        Set<String> unClean = getRegions();
        Set<String> result = new HashSet<String>();
        for ( String region : unClean ) {
            String newR = region.toLowerCase().trim();
            newR = newR.replace( "_", "" );
            result.add( newR );
        }
        return result;
    }

    /**
     * returns all Allen regions
     */
    public Set<String> getRegions() {
        // assumes all names have acronyms
        return new HashSet<String>( acronyms.values() );
    }

    public Set<String> getBAMSMappedRegions() {
        return BAMSMap.getSeenValues();
    }

    public Set<String> getAllenMappedRegions() {
        return BAMSMap.keySet();
    }

    public static void main( String args[] ) throws Exception {
        StructureCatalogLoader test = new StructureCatalogLoader();

        System.out.println( test.getParent( "Nucleus x" ) );
        System.out.println( test.getParents( "Nucleus x" ) );

        System.out.println( test.getFullName( "PD" ) );
        System.out.println( test.getFullName( "x" ) );
        System.out.println( test.getFullName( "y" ) );

        System.out.println( test.getParent( "Posterodorsal preoptic nucleus" ) );
        System.out.println( test.getParent( "Nucleus x" ) );
        log.info( test.getParent( "Nucleus y" ) );

        log.info( test.getChildren( "Cerebral nuclei" ) );

        // test.doStats();

        log.info( test.getNomenclatureMatrix() );

        log.info( test.getMotorNuclei() );
        // System.out.println( "BAMS leafs:" + bams.getLeafAmount() );
        // Set<String> bamsleafs = bams.getRegions();
        // for ( String region : test.getLeafs() ) {
        // if ( test.getBAMSMappedRegions( region ) != null ) {
        // bamsleafs.removeAll( test.getBAMSMappedRegions( region ) );
        // }
        // }
        // System.out.println( "BAMS leafs minus Allen leafs:" + bamsleafs.size() );
        // System.out.println(bamsleafs);

        // bams.getRegionAmount() +","+ );

    }

    public Set<String> getMotorNuclei() {
        return motorNuclei;
    }
}
