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
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

import javax.vecmath.Point3d;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import ubic.BAMSandAllen.SetupParameters;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;
import au.com.bytecode.opencsv.CSVReader;

public class AllenAtlasAnnotationLoader {
    private static Log log = LogFactory.getLog( AllenAtlasAnnotationLoader.class.getName() );
    Map<String, Set<Point3d>> atlasMap;
    BrainStructuresCSVLoader apiLoader;

    // Like the other volume files available through this API, the annotated atlases consist of a series of coordinates,
    // each with an associated value. In this case the value at each coordinate point is a brain structure identifier.

    // The header section of the annotated atlas volumes lists the resolution in the comments section:
    //
    // Comment:Atlas_200
    // Dimensions:67, 41, 58
    //
    // The data value associated with each coordinate is the Informatics ID of the brain structure at that point in the
    // volume. See the Brain Structure Ontology section for details on the structure IDs & names.
    public AllenAtlasAnnotationLoader() throws Exception {
        this( SetupParameters.config.getString( "abams.allenAtlasFile" ) );
    }

    public AllenAtlasAnnotationLoader( String filename ) throws Exception {
        StructureCatalogLoader mainLoader = new StructureCatalogLoader();
        apiLoader = new BrainStructuresCSVLoader();

        CSVReader reader = new CSVReader( new FileReader( filename ) );
        // List<String[]> fileContents = reader.readAll();
        String[] line;
        atlasMap = new HashMap<String, Set<Point3d>>();
        while ( null != ( line = reader.readNext() ) ) {
            if ( line[0].startsWith( "Comment" ) || line[0].startsWith( "Dimensions" ) ) {
                continue;
            }
            int x = Integer.parseInt( line[0] );
            int y = Integer.parseInt( line[1] );
            int z = Integer.parseInt( line[2] );
            String id = line[3];
            Point3d point = new Point3d( x, y, z );
            addPoint( id, point );
            // get all parents and add this point to them?
            String name = apiLoader.getName( id );
            Set<String> parents = mainLoader.getParents( name );
            // log.info( "adding parents of " + name + ":" + parents );
            for ( String parent : parents ) {
                String parentID = apiLoader.getID( parent );
                addPoint( parentID, point );
            }
        }
    }

    public Set<Point3d> getVoxels( String regionName ) {
        String parentID = apiLoader.getID( regionName );
        return atlasMap.get( parentID );
    }

    public int getTotalVolume( Collection<String> regions ) {
        Set<Point3d> total = new HashSet<Point3d>();
        for ( String id : atlasMap.keySet() ) {
            String name = apiLoader.getName( id );
            if ( regions == null || regions.contains( name ) ) {
                total.addAll( atlasMap.get( id ) );
            }
        }
        return total.size();
    }

    private void addPoint( String id, Point3d point ) {
        Set<Point3d> points = atlasMap.get( id );
        if ( points == null ) {
            points = new HashSet<Point3d>();
            atlasMap.put( id, points );
        }
        points.add( point );
    }

    public void addPoints( String id, Set<Point3d> points ) {
        for ( Point3d point : points ) {
            addPoint( id, point );
        }
    }

    public Map<String, Point3d> getDimensions() {
        Map<String, Point3d> result = new HashMap<String, Point3d>();
        for ( String id : atlasMap.keySet() ) {
            double maxx, minx, miny, minz, maxy, maxz;
            maxx = maxy = maxz = Double.MIN_VALUE;
            minx = miny = minz = Double.MAX_VALUE;
            Set<Point3d> points = atlasMap.get( id );
            // find min and max of x,y,z
            for ( Point3d p : points ) {
                maxx = Math.max( maxx, p.x );
                maxy = Math.max( maxy, p.y );
                maxz = Math.max( maxz, p.z );
                minx = Math.min( minx, p.x );
                miny = Math.min( miny, p.y );
                minz = Math.min( minz, p.z );
            }
            // compute dimensions
            String name = apiLoader.getName( id );

            Point3d dims = new Point3d( maxx - minx, maxy - miny, maxz - minz );
            // log.info(name+"->"+ dims);
            result.put( name, dims );
        }
        return result;
    }

    public Map<String, Point3d> getCenters() {
        Map<String, Point3d> result = new HashMap<String, Point3d>();
        for ( String id : atlasMap.keySet() ) {
            Set<Point3d> points = atlasMap.get( id );
            Point3d centre = new Point3d( 0, 0, 0 );
            for ( Point3d point : points ) {
                centre.add( point );
            }
            // log.info( points.size() + " id:" + id );
            centre.scale( 1d / points.size() );
            String name = apiLoader.getName( id );
            // log.info( name + " at " + centre.toString() );
            result.put( name, centre );
        }
        return result;
    }

    public DoubleMatrix<String, String> getDimensionsMatrix() {
        Map<String, Point3d> dims = getDimensions();
        return getPoints2Matrix( dims );
    }

    /**
     * Computes a matrix of brain region volumes by using the voxel counts in the atlas.
     * 
     * @return
     */
    public DoubleMatrix<String, String> getVolumeMatrix() {
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( 1, apiLoader.getNames().size() );

        result.setColumnNames( new LinkedList<String>( apiLoader.getNames() ) );
        result.addRowName( "volume" );

        for ( String region : atlasMap.keySet() ) {
            String regionName = apiLoader.getName( region );
            double voxels = atlasMap.get( region ).size();
            result.setByKeys( "volume", regionName, voxels );
        }
        return result;
    }

    /**
     * return a matrix where the column names are the allen regions and the rows are the coordinates of the centers of
     * each region
     */
    public DoubleMatrix<String, String> getCenterMatrix() {
        Map<String, Point3d> centers = getCenters();
        return getPoints2Matrix( centers );
    }

    public DoubleMatrix<String, String> getPoints2Matrix( Map<String, Point3d> points ) {
        DoubleMatrix<String, String> result = new DenseDoubleMatrix<String, String>( 3, points.size() );
        result.setColumnNames( new LinkedList<String>( points.keySet() ) );
        result.addRowName( "x" );
        result.addRowName( "y" );
        result.addRowName( "z" );
        for ( String region : points.keySet() ) {
            Point3d center = points.get( region );
            result.setByKeys( "x", region, center.x );
            result.setByKeys( "y", region, center.y );
            result.setByKeys( "z", region, center.z );
        }
        // log.info( result );
        // log.info( "Rows:" + result.rows() );
        // log.info( "Cols:" + result.columns() );
        return result;
    }

    /**
     * @param args
     */
    public static void main( String[] args ) throws Exception {
        // TODO Auto-generated method stub
        AllenAtlasAnnotationLoader atlas = new AllenAtlasAnnotationLoader();
        atlas.getCenters();
        atlas.getCenterMatrix();
        log.info( atlas.getTotalVolume( null ) );

        log.info( atlas.getVoxels( "Interstitial nucleus of Cajal" ) );
        log.info( atlas.getVoxels( "Nucleus of Darkschewitsch" ) );
    }
}
