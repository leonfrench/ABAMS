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
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import ubic.BAMSandAllen.SetupParameters;
import ubic.BAMSandAllen.Util;
import ubic.basecode.dataStructure.CountingMap;
import ubic.basecode.dataStructure.matrix.DenseDoubleMatrix;
import ubic.basecode.dataStructure.matrix.DoubleMatrix;

public class BAMSXMLConnectivityLoader {
    protected static Log log = LogFactory.getLog( BAMSXMLConnectivityLoader.class );
    boolean generateMatrix = false;
    Set<String> regionNames;
    Set<String> uniqueConnections;
    Set<String> uniqueNotPresentConnections;

    DoubleMatrix<String, String> outgoingMatrix;
    DoubleMatrix<String, String> notPresentMatrix;

    CountingMap<String> strengths;
    int processConnectionCalls = 0;

    public BAMSXMLConnectivityLoader() throws Exception {
        generateMatrix = false;
        BAMSDataLoader swanson98 = new BAMSDataLoader();
        strengths = new CountingMap<String>();
        uniqueConnections = new HashSet<String>();
        uniqueNotPresentConnections = new HashSet<String>();
        regionNames = new HashSet<String>();

        DocumentBuilderFactory docBuilderFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder docBuilder = docBuilderFactory.newDocumentBuilder();
        String filename = SetupParameters.config.getString( "whitetext.lexicon.output" )
                + "bams-connectome-october-2011.xml";
        Document doc = docBuilder.parse( new File( filename ) );

        // normalize text representation
        doc.getDocumentElement().normalize();
        System.out.println( "Root element of the doc is " + doc.getDocumentElement().getNodeName() );

        NodeList gray_matter = doc.getElementsByTagName( "gray_matter_region" );
        int totalRegions = gray_matter.getLength();
        System.out.println( "Total no of regions : " + totalRegions );
        // System.exit( 1 );

        processXML( gray_matter );

        log.info( "Names: " + regionNames.size() );
        Set<String> swanson98Regions = swanson98.getRegions();
        log.info( "Intersection of swanson 98 and this xml file regions: "
                + Util.intersectSize( swanson98Regions, regionNames ) );

        generateMatrix = true;
        // now that I know how many regions I can make the connection matrix
        outgoingMatrix = new DenseDoubleMatrix<String, String>( regionNames.size(), regionNames.size() );
        notPresentMatrix = new DenseDoubleMatrix<String, String>( regionNames.size(), regionNames.size() );
        outgoingMatrix.setRowNames( new LinkedList<String>( regionNames ) );
        outgoingMatrix.setColumnNames( new LinkedList<String>( regionNames ) );
        notPresentMatrix.setRowNames( new LinkedList<String>( regionNames ) );
        notPresentMatrix.setColumnNames( new LinkedList<String>( regionNames ) );

        // reset before reloading connections
        strengths = new CountingMap<String>();

        processXML( gray_matter );

    }

    private void processXML( NodeList gray_matter ) {
        for ( int s = 0; s < gray_matter.getLength(); s++ ) {

            Node firstRegionNode = gray_matter.item( s );
            if ( firstRegionNode.getNodeType() == Node.ELEMENT_NODE ) {

                Element firstPersonElement = ( Element ) firstRegionNode;

                // -------
                NodeList firstNameList = firstPersonElement.getElementsByTagName( "name" );
                Element firstNameElement = ( Element ) firstNameList.item( 0 );

                NodeList textFNList = firstNameElement.getChildNodes();
                String regionName = ( ( Node ) textFNList.item( 0 ) ).getNodeValue().trim();
                // System.out.println( "Name : " + regionName );

                regionNames.add( regionName );

                // targets
                NodeList targets = firstPersonElement.getElementsByTagName( "targets" );
                NodeList sources = firstPersonElement.getElementsByTagName( "sources" );

                getConnectionsFromBlock( regionName, targets );
                getConnectionsFromBlock( regionName, sources );

            }// end of if clause

        }// end of for loop with s var
    }

    private void getConnectionsFromBlock( String regionName, NodeList targets ) {
        Element tag = ( Element ) targets.item( 0 );
        if ( tag != null ) {
            NodeList subNodes = tag.getChildNodes();
            String strength = null;
            String connectedName = null;
            boolean outgoing = false;
            String abbreviation = null;
            for ( int ss = 0; ss < subNodes.getLength(); ss++ ) {
                Node connectionSubNode = subNodes.item( ss );

                if ( connectionSubNode.getNodeType() == Node.ELEMENT_NODE ) {
                    String tagValue = connectionSubNode.getTextContent().trim();
                    String tagName = connectionSubNode.getNodeName();
                    // log.info( tagName + " = " + tagValue );

                    if ( tagName.startsWith( "source" ) ) outgoing = false;
                    if ( tagName.startsWith( "target" ) ) outgoing = true;

                    if ( tagName.equals( "targetname" ) || tagName.equals( "sourcename" ) ) {
                        connectedName = tagValue;
                    }
                    if ( tagName.equals( "targetabbreviation" ) || tagName.equals( "sourceabbreviation" ) ) {
                        abbreviation = tagValue;
                    }
                    if ( tagName.equals( "strength" ) ) {
                        strength = tagValue;
                        processConnection( regionName, outgoing, strength, connectedName, abbreviation );

                        strength = null;
                        connectedName = null;
                        abbreviation = null;
                    }
                }
                // targetSubNode.getChildNodes();
                // String targetSubNodeName = ( ( Node ) textFNList.item( 0 ) ).getNodeValue().trim();
            }
        } // end if to check for any connections
    }

    public void processConnection( String regionName, boolean outgoing, String strength, String targetname,
            String targetabbreviation ) {
        processConnectionCalls++;
        strengths.increment( strength );
        regionNames.add( targetname );

        if ( !outgoing ) {
            String temp = regionName;
            regionName = targetname;
            targetname = temp;
        }

        // log.info( regionName + targetname );
        if ( strength.equals( "not present" ) || strength.equals( "fibers of passage" ) ) {
            uniqueNotPresentConnections.add( regionName + targetname );
            if ( generateMatrix ) notPresentMatrix.setByKeys( regionName, targetname, 1d );
        } else {
            uniqueConnections.add( regionName + targetname );
            if ( generateMatrix ) outgoingMatrix.setByKeys( regionName, targetname, 1d );
        }

    }

    public void printStats() {
        // log.info( "Unique connections/2:" + uniqueConnections.size() / 2 );
        log.info( "Unique connections:" + uniqueConnections.size() );
        log.info( "Unique not present connections:" + uniqueNotPresentConnections.size() );
        System.out.println( strengths.toString() );
        System.out.println( "Number of strenghts:" + strengths.summation() );
        log.info( "Calls to processconnection:" + processConnectionCalls );
        log.info( "Both present and not:" + Util.intersectSize( uniqueConnections, uniqueNotPresentConnections ) );
    }

    public static void main( String args[] ) throws Exception {
        BAMSXMLConnectivityLoader loader = new BAMSXMLConnectivityLoader();
        loader.printStats();

    }

    public DoubleMatrix<String, String> getNotPresentMatrix() {
        return notPresentMatrix;
    }

    public DoubleMatrix<String, String> getOutgoingMatrix() {
        return outgoingMatrix;
    }

}
