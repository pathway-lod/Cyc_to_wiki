from scripts.data_structure.wiki_data_structure import (
    DataNode, Interaction, Pathway, HAlign, VAlign, BorderStyle, ShapeType, LineStyle,
    ConnectorType, ArrowHeadType, AnchorShapeType, Point, Anchor, Graphics, Citation
)


class GPMLWriter:
    def __init__(self):
        pass

    def escape_xml(self, text):
        """
        Escape special XML characters in text

        Args:
            text: String to escape

        Returns:
            str: Escaped string safe for XML attributes
        """
        if text is None:
            return ""
        if not isinstance(text, str):
            text = str(text)


        text = text.replace('&', '&amp;')
        text = text.replace('<', '&lt;')
        text = text.replace('>', '&gt;')
        text = text.replace('"', '&quot;')
        text = text.replace("'", '&apos;')

        return text

    def ensure_interaction_graphics(self, interaction):
        """
        Ensure an interaction has valid graphics with all required attributes

        Args:
            interaction: Interaction object

        Returns:
            Graphics object with required defaults
        """
        if not hasattr(interaction, 'graphics') or interaction.graphics is None:
            # Create default graphics
            interaction.graphics = Graphics(
                lineColor='000000',
                lineStyle=LineStyle.SOLID,
                lineWidth=1.0,
                connectorType=ConnectorType.STRAIGHT
            )
        else:
            if interaction.graphics.lineColor is None:
                interaction.graphics.lineColor = '000000'
            if interaction.graphics.lineStyle is None:
                interaction.graphics.lineStyle = LineStyle.SOLID
            if interaction.graphics.lineWidth is None:
                interaction.graphics.lineWidth = 1.0
            if interaction.graphics.connectorType is None:
                interaction.graphics.connectorType = ConnectorType.STRAIGHT

        return interaction.graphics


    def write_datanode(self, datanode: DataNode) -> str:
        """
        Converts a DataNode object into a GPML DataNode XML string.

        Args:
            datanode: The DataNode object to convert.

        Returns:
            str: The GPML XML string for the DataNode.
        """
        element_id = self.escape_xml(datanode.elementId)
        text_label = self.escape_xml(datanode.textLabel)
        node_type = datanode.type.value

        gpml_output = f'    <DataNode elementId="{element_id}" textLabel="{text_label}" type="{node_type}"'

        if hasattr(datanode, 'groupRef') and datanode.groupRef:
            gpml_output += f' groupRef="{self.escape_xml(datanode.groupRef)}"'
        if hasattr(datanode, 'aliasRef') and datanode.aliasRef:
            gpml_output += f' aliasRef="{self.escape_xml(datanode.aliasRef)}"'

        gpml_output += '>\n'

        # Handle Xref - prioritize single xref over xrefs list
        if hasattr(datanode, 'xref') and datanode.xref:
            gpml_output += f'      <Xref identifier="{self.escape_xml(datanode.xref.identifier)}" dataSource="{self.escape_xml(datanode.xref.dataSource)}" />\n'
        elif hasattr(datanode, 'xrefs') and datanode.xrefs and len(datanode.xrefs) > 0:
            # Use the first xref as primary
            xref = datanode.xrefs[0]
            gpml_output += f'      <Xref identifier="{self.escape_xml(xref.identifier)}" dataSource="{self.escape_xml(xref.dataSource)}" />\n'

        # States (if present)
        if hasattr(datanode, 'states') and datanode.states:
            gpml_output += '      <States>\n'
            for state in datanode.states:
                gpml_output += self.write_state(state)
            gpml_output += '      </States>\n'

        # Graphics section - ensure all required attributes are present
        gpml_output += '      <Graphics'

        # Required position and dimension attributes
        gpml_output += f' centerX="{datanode.graphics.centerX}"'
        gpml_output += f' centerY="{datanode.graphics.centerY}"'
        gpml_output += f' width="{datanode.graphics.width}"'
        gpml_output += f' height="{datanode.graphics.height}"'

        # Text color
        if datanode.graphics.textColor is not None:
            gpml_output += f' textColor="{datanode.graphics.textColor}"'
        else:
            gpml_output += ' textColor="000000"'

        # Font attributes
        # FontName
        font_name = datanode.graphics.fontName if datanode.graphics.fontName is not None else "Arial"
        gpml_output += f' fontName="{self.escape_xml(font_name)}"'

        # FontWeight
        if datanode.graphics.fontWeight is not None:
            gpml_output += ' fontWeight="Bold"' if datanode.graphics.fontWeight else ' fontWeight="Normal"'
        else:
            gpml_output += ' fontWeight="Normal"'

        # FontStyle
        if datanode.graphics.fontStyle is not None:
            gpml_output += ' fontStyle="Italic"' if datanode.graphics.fontStyle else ' fontStyle="Normal"'
        else:
            gpml_output += ' fontStyle="Normal"'

        # FontDecoration
        if datanode.graphics.fontDecoration is not None:
            gpml_output += ' fontDecoration="Underline"' if datanode.graphics.fontDecoration else ' fontDecoration="Normal"'
        else:
            gpml_output += ' fontDecoration="Normal"'

        # FontStrikethru
        if datanode.graphics.fontStrikethru is not None:
            gpml_output += ' fontStrikethru="Strikethru"' if datanode.graphics.fontStrikethru else ' fontStrikethru="Normal"'
        else:
            gpml_output += ' fontStrikethru="Normal"'

        # FontSize
        font_size = int(datanode.graphics.fontSize) if datanode.graphics.fontSize is not None else 10
        gpml_output += f' fontSize="{font_size}"'

        # Alignment
        if datanode.graphics.hAlign is not None:
            gpml_output += f' hAlign="{datanode.graphics.hAlign.value}"'
        else:
            gpml_output += ' hAlign="Center"'

        if datanode.graphics.vAlign is not None:
            gpml_output += f' vAlign="{datanode.graphics.vAlign.value}"'
        else:
            gpml_output += ' vAlign="Middle"'

        # Shape style attributes
        # Border color
        border_color = datanode.graphics.borderColor if datanode.graphics.borderColor is not None else "000000"
        gpml_output += f' borderColor="{border_color}"'

        # Border style
        if datanode.graphics.borderStyle is not None:
            gpml_output += f' borderStyle="{datanode.graphics.borderStyle.value}"'
        else:
            gpml_output += ' borderStyle="Solid"'

        # Border width
        border_width = datanode.graphics.borderWidth if datanode.graphics.borderWidth is not None else 1.0
        gpml_output += f' borderWidth="{border_width}"'

        # Fill color
        fill_color = datanode.graphics.fillColor if datanode.graphics.fillColor is not None else "FFFFFF"
        gpml_output += f' fillColor="{fill_color}"'

        # Shape type
        if datanode.graphics.shapeType is not None:
            gpml_output += f' shapeType="{datanode.graphics.shapeType.value}"'
        else:
            gpml_output += ' shapeType="Rectangle"'

        # Optional attributes
        if datanode.graphics.zOrder is not None:
            gpml_output += f' zOrder="{int(datanode.graphics.zOrder)}"'

        if hasattr(datanode.graphics, 'rotation') and datanode.graphics.rotation is not None:
            gpml_output += f' rotation="{datanode.graphics.rotation}"'

        gpml_output += ' />\n'

        # Comments
        if hasattr(datanode, 'comments') and datanode.comments:
            for comment in datanode.comments:
                if isinstance(comment, dict):
                    source = comment.get("source", "")
                    text = comment.get("text", comment.get("value", ""))
                else:
                    source = getattr(comment, 'source', '')
                    text = getattr(comment, 'value', '')

                if source:
                    gpml_output += f'      <Comment source="{self.escape_xml(source)}">{self.escape_xml(text)}</Comment>\n'
                else:
                    gpml_output += f'      <Comment>{self.escape_xml(text)}</Comment>\n'

        # Properties
        if hasattr(datanode, 'properties') and datanode.properties:
            for prop in datanode.properties:
                if isinstance(prop, dict):
                    key = prop.get("key", "")
                    value = prop.get("value", "")
                else:
                    key = getattr(prop, 'key', '')
                    value = getattr(prop, 'value', '')

                gpml_output += f'      <Property key="{self.escape_xml(key)}" value="{self.escape_xml(value)}" />\n'

        # AnnotationRefs
        if hasattr(datanode, 'annotationRefs') and datanode.annotationRefs:
            for annotation_ref in datanode.annotationRefs:
                ref_id = getattr(annotation_ref, 'elementRef', '')
                gpml_output += f'      <AnnotationRef elementRef="{self.escape_xml(ref_id)}" />\n'

        # CitationRefs
        if hasattr(datanode, 'citationRefs') and datanode.citationRefs:
            for citation_ref in datanode.citationRefs:
                ref_id = getattr(citation_ref, 'elementRef', '')
                gpml_output += f'      <CitationRef elementRef="{self.escape_xml(ref_id)}" />\n'

        # EvidenceRefs
        if hasattr(datanode, 'evidenceRefs') and datanode.evidenceRefs:
            for evidence_ref in datanode.evidenceRefs:
                ref_id = getattr(evidence_ref, 'elementRef', '')
                gpml_output += f'      <EvidenceRef elementRef="{self.escape_xml(ref_id)}" />\n'

        gpml_output += '    </DataNode>\n'
        return gpml_output

    def write_group(self, group) -> str:
        """
        Write a Group element to GPML format.
        Groups are visual containers for complex components.

        Args:
            group: Group object representing a complex or other grouping

        Returns:
            str: GPML XML string for the Group
        """
        element_id = self.escape_xml(group.elementId)
        gpml_output = f'    <Group elementId="{element_id}"'

        # Optional textLabel
        if hasattr(group, 'textLabel') and group.textLabel:
            gpml_output += f' textLabel="{self.escape_xml(group.textLabel)}"'

        # Group type
        if hasattr(group, 'type') and group.type:
            gpml_output += f' type="{group.type.value}"'

        # GroupRef for nested groups (usually not needed for complexes)
        if hasattr(group, 'groupRef') and group.groupRef:
            gpml_output += f' groupRef="{self.escape_xml(group.groupRef)}"'

        gpml_output += '>\n'

        # Xref if present
        if hasattr(group, 'xref') and group.xref:
            gpml_output += f'      <Xref identifier="{self.escape_xml(group.xref.identifier)}" dataSource="{self.escape_xml(group.xref.dataSource)}" />\n'

        # Graphics section - REQUIRED for Groups
        if hasattr(group, 'graphics') and group.graphics:
            gpml_output += '      <Graphics'

            # Required position and dimension attributes
            gpml_output += f' centerX="{group.graphics.centerX}"'
            gpml_output += f' centerY="{group.graphics.centerY}"'
            gpml_output += f' width="{group.graphics.width}"'
            gpml_output += f' height="{group.graphics.height}"'

            # Text color (for label)
            if group.graphics.textColor is not None:
                gpml_output += f' textColor="{group.graphics.textColor}"'
            else:
                gpml_output += ' textColor="000000"'

            # Font attributes
            if group.graphics.fontName:
                gpml_output += f' fontName="{self.escape_xml(group.graphics.fontName)}"'

            if group.graphics.fontSize is not None:
                gpml_output += f' fontSize="{int(group.graphics.fontSize)}"'

            # Border and fill - important for visual representation
            if group.graphics.borderColor:
                gpml_output += f' borderColor="{group.graphics.borderColor}"'
            else:
                gpml_output += ' borderColor="999999"'  # Default gray border

            if group.graphics.borderStyle:
                gpml_output += f' borderStyle="{group.graphics.borderStyle.value}"'
            else:
                gpml_output += ' borderStyle="Dashed"'  # Default dashed for groups

            if group.graphics.borderWidth is not None:
                gpml_output += f' borderWidth="{group.graphics.borderWidth}"'
            else:
                gpml_output += ' borderWidth="2"'

            # Transparent or light fill for groups
            if group.graphics.fillColor:
                gpml_output += f' fillColor="{group.graphics.fillColor}"'
            else:
                gpml_output += ' fillColor="F0F0F0"'  # Light gray default

            if group.graphics.shapeType:
                gpml_output += f' shapeType="{group.graphics.shapeType.value}"'
            else:
                gpml_output += ' shapeType="Rectangle"'

            if group.graphics.zOrder is not None:
                gpml_output += f' zOrder="{int(group.graphics.zOrder)}"'
            else:
                gpml_output += ' zOrder="0"'  # Groups typically behind other elements

            gpml_output += ' />\n'

        # Comments
        if hasattr(group, 'comments') and group.comments:
            for comment in group.comments:
                if isinstance(comment, dict):
                    source = comment.get("source", "")
                    text = comment.get("text", comment.get("value", ""))
                else:
                    source = getattr(comment, 'source', '')
                    text = getattr(comment, 'value', '')

                if source:
                    gpml_output += f'      <Comment source="{self.escape_xml(source)}">{self.escape_xml(text)}</Comment>\n'
                else:
                    gpml_output += f'      <Comment>{self.escape_xml(text)}</Comment>\n'

        # Properties
        if hasattr(group, 'properties') and group.properties:
            for prop in group.properties:
                if isinstance(prop, dict):
                    key = prop.get("key", "")
                    value = prop.get("value", "")
                else:
                    key = getattr(prop, 'key', '')
                    value = getattr(prop, 'value', '')

                gpml_output += f'      <Property key="{self.escape_xml(key)}" value="{self.escape_xml(value)}" />\n'

        # CitationRefs
        if hasattr(group, 'citationRefs') and group.citationRefs:
            for citation_ref in group.citationRefs:
                ref_id = getattr(citation_ref, 'elementRef', '')
                gpml_output += f'      <CitationRef elementRef="{self.escape_xml(ref_id)}" />\n'

        gpml_output += '    </Group>\n'
        return gpml_output

    def write_state(self, state) -> str:
        """
        Write a State element for a DataNode.

        Args:
            state: State object

        Returns:
            str: GPML XML string for the State
        """
        element_id = self.escape_xml(state.elementId)
        text_label = self.escape_xml(state.textLabel)

        gpml_output = f'        <State elementId="{element_id}" textLabel="{text_label}"'

        if hasattr(state, 'type') and state.type:
            gpml_output += f' type="{state.type.value}"'

        gpml_output += '>\n'

        # Xref if present
        if hasattr(state, 'xref') and state.xref:
            gpml_output += f'          <Xref identifier="{self.escape_xml(state.xref.identifier)}" dataSource="{self.escape_xml(state.xref.dataSource)}" />\n'

        # Graphics - required for states
        if hasattr(state, 'graphics') and state.graphics:
            gpml_output += '          <Graphics'

            # Required attributes for state graphics
            gpml_output += f' relX="{state.graphics.relX}"'
            gpml_output += f' relY="{state.graphics.relY}"'
            gpml_output += f' width="{state.graphics.width}"'
            gpml_output += f' height="{state.graphics.height}"'

            # Optional font attributes
            if state.graphics.textColor:
                gpml_output += f' textColor="{state.graphics.textColor}"'

            font_name = state.graphics.fontName if state.graphics.fontName else "Arial"
            gpml_output += f' fontName="{self.escape_xml(font_name)}"'

            if state.graphics.fontWeight is not None:
                gpml_output += ' fontWeight="Bold"' if state.graphics.fontWeight else ' fontWeight="Normal"'

            if state.graphics.fontStyle is not None:
                gpml_output += ' fontStyle="Italic"' if state.graphics.fontStyle else ' fontStyle="Normal"'

            font_size = int(state.graphics.fontSize) if state.graphics.fontSize else 10
            gpml_output += f' fontSize="{font_size}"'

            # Shape attributes
            if state.graphics.borderColor:
                gpml_output += f' borderColor="{state.graphics.borderColor}"'

            if state.graphics.borderStyle:
                gpml_output += f' borderStyle="{state.graphics.borderStyle.value}"'

            if state.graphics.borderWidth is not None:
                gpml_output += f' borderWidth="{state.graphics.borderWidth}"'

            if state.graphics.fillColor:
                gpml_output += f' fillColor="{state.graphics.fillColor}"'

            if state.graphics.shapeType:
                gpml_output += f' shapeType="{state.graphics.shapeType.value}"'

            gpml_output += ' />\n'

        # Comments
        if hasattr(state, 'comments') and state.comments:
            for comment in state.comments:
                source = getattr(comment, 'source', '')
                text = getattr(comment, 'value', '')
                if source:
                    gpml_output += f'          <Comment source="{self.escape_xml(source)}">{self.escape_xml(text)}</Comment>\n'
                else:
                    gpml_output += f'          <Comment>{self.escape_xml(text)}</Comment>\n'

        # Properties
        if hasattr(state, 'properties') and state.properties:
            for prop in state.properties:
                key = getattr(prop, 'key', '')
                value = getattr(prop, 'value', '')
                gpml_output += f'          <Property key="{self.escape_xml(key)}" value="{self.escape_xml(value)}" />\n'

        # CitationRefs
        if hasattr(state, 'citationRefs') and state.citationRefs:
            for citation_ref in state.citationRefs:
                ref_id = getattr(citation_ref, 'elementRef', '')
                gpml_output += f'          <CitationRef elementRef="{self.escape_xml(ref_id)}" />\n'

        gpml_output += '        </State>\n'
        return gpml_output

    def write_interaction(self, interaction: Interaction) -> str:
        """
        Converts an Interaction object into a GPML Interaction XML string.
        IMPORTANT: Graphics element is REQUIRED for Interactions in GPML 2021!

        Args:
            interaction: The Interaction object to convert.

        Returns:
            str: The GPML XML string for the Interaction.
        """
        # Ensure graphics exists with required attributes
        self.ensure_interaction_graphics(interaction)

        element_id = self.escape_xml(interaction.elementId)
        gpml_output = f'    <Interaction elementId="{element_id}"'
        if hasattr(interaction, 'groupRef') and interaction.groupRef:
            gpml_output += f' groupRef="{self.escape_xml(interaction.groupRef)}"'

        gpml_output += '>\n'

        # Write Xref
        if interaction.xref:
            gpml_output += f'      <Xref identifier="{self.escape_xml(interaction.xref.identifier)}" dataSource="{self.escape_xml(interaction.xref.dataSource)}" />\n'

        # Write Waypoints
        gpml_output += '      <Waypoints>\n'

        for point in interaction.waypoints:
            gpml_output += '        <Point'
            gpml_output += f' elementId="{self.escape_xml(point.elementId)}"'

            # ArrowHead is optional, only add if not default
            if hasattr(point, 'arrowHead') and point.arrowHead != ArrowHeadType.UNDIRECTED:
                gpml_output += f' arrowHead="{point.arrowHead.value}"'

            # Coordinates
            gpml_output += f' x="{point.x}" y="{point.y}"'

            # Optional attributes
            if hasattr(point, 'elementRef') and point.elementRef:
                gpml_output += f' elementRef="{self.escape_xml(point.elementRef)}"'
            if hasattr(point, 'relX') and point.relX is not None:
                gpml_output += f' relX="{point.relX}"'
            if hasattr(point, 'relY') and point.relY is not None:
                gpml_output += f' relY="{point.relY}"'

            gpml_output += ' />\n'

        # Write Anchors if present
        if hasattr(interaction, 'anchors') and interaction.anchors:
            for anchor in interaction.anchors:
                gpml_output += '        <Anchor'
                gpml_output += f' elementId="{self.escape_xml(anchor.elementId)}"'
                gpml_output += f' position="{anchor.position}"'

                if hasattr(anchor, 'shapeType') and anchor.shapeType != AnchorShapeType.SQUARE:
                    gpml_output += f' shapeType="{anchor.shapeType.value}"'

                gpml_output += ' />\n'

        gpml_output += '      </Waypoints>\n'

        # Write Graphics
        gpml_output += '      <Graphics'
        gpml_output += f' lineColor="{interaction.graphics.lineColor}"'
        gpml_output += f' lineStyle="{interaction.graphics.lineStyle.value}"'
        gpml_output += f' lineWidth="{interaction.graphics.lineWidth}"'
        gpml_output += f' connectorType="{interaction.graphics.connectorType.value}"'

        if hasattr(interaction.graphics, 'zOrder') and interaction.graphics.zOrder is not None:
            gpml_output += f' zOrder="{int(interaction.graphics.zOrder)}"'

        gpml_output += ' />\n'

        # Comments
        if hasattr(interaction, 'comments') and interaction.comments:
            for comment in interaction.comments:
                if isinstance(comment, dict):
                    source = comment.get("source", "")
                    text = comment.get("text", comment.get("value", ""))
                else:
                    source = getattr(comment, 'source', '')
                    text = getattr(comment, 'value', '')

                if source:
                    gpml_output += f'      <Comment source="{self.escape_xml(source)}">{self.escape_xml(text)}</Comment>\n'
                else:
                    gpml_output += f'      <Comment>{self.escape_xml(text)}</Comment>\n'

        # Properties
        if hasattr(interaction, 'properties') and interaction.properties:
            for prop in interaction.properties:
                if isinstance(prop, dict):
                    key = prop.get("key", "")
                    value = prop.get("value", "")
                else:
                    key = getattr(prop, 'key', '')
                    value = getattr(prop, 'value', '')

                gpml_output += f'      <Property key="{self.escape_xml(key)}" value="{self.escape_xml(value)}" />\n'

        # AnnotationRefs
        if hasattr(interaction, 'annotationRefs') and interaction.annotationRefs:
            for annotation_ref in interaction.annotationRefs:
                ref_id = getattr(annotation_ref, 'elementRef', '')
                gpml_output += f'      <AnnotationRef elementRef="{self.escape_xml(ref_id)}" />\n'

        # CitationRefs
        if hasattr(interaction, 'citationRefs') and interaction.citationRefs:
            for citation_ref in interaction.citationRefs:
                ref_id = getattr(citation_ref, 'elementRef', '')
                gpml_output += f'      <CitationRef elementRef="{self.escape_xml(ref_id)}" />\n'

        # EvidenceRefs
        if hasattr(interaction, 'evidenceRefs') and interaction.evidenceRefs:
            for evidence_ref in interaction.evidenceRefs:
                ref_id = getattr(evidence_ref, 'elementRef', '')
                gpml_output += f'      <EvidenceRef elementRef="{self.escape_xml(ref_id)}" />\n'

        gpml_output += '    </Interaction>\n'
        return gpml_output

    def write_pathway(self, pathway: Pathway) -> str:
        """
        Write a complete pathway to GPML format with proper Groups support

        Args:
            pathway: Pathway object containing all elements

        Returns:
            str: Complete GPML XML string
        """
        # Required XML declaration
        gpml_output = '<?xml version="1.0" encoding="UTF-8"?>\n'

        # Root Pathway element
        gpml_output += '<Pathway xmlns="http://pathvisio.org/GPML/2021"'

        # Title attribute
        gpml_output += f' title="{self.escape_xml(pathway.title)}"'

        # Optional attributes
        if pathway.organism:
            gpml_output += f' organism="{self.escape_xml(pathway.organism)}"'
        if pathway.source:
            gpml_output += f' source="{self.escape_xml(pathway.source)}"'
        if pathway.version:
            gpml_output += f' version="{self.escape_xml(pathway.version)}"'
        if pathway.license:
            gpml_output += f' license="{self.escape_xml(pathway.license)}"'

        gpml_output += '>\n'

        # Pathway Xref
        if pathway.xref:
            gpml_output += f'  <Xref identifier="{self.escape_xml(pathway.xref.identifier)}" dataSource="{self.escape_xml(pathway.xref.dataSource)}" />\n'

        # Description
        if pathway.description:
            gpml_output += f'  <Description>{self.escape_xml(pathway.description)}</Description>\n'

        # Authors
        if pathway.authors:
            gpml_output += '  <Authors>\n'
            for author in pathway.authors:
                gpml_output += f'    <Author name="{self.escape_xml(author.name)}"'

                if hasattr(author, 'username') and author.username:
                    gpml_output += f' username="{self.escape_xml(author.username)}"'
                if hasattr(author, 'order') and author.order is not None:
                    gpml_output += f' order="{author.order}"'

                if hasattr(author, 'xref') and author.xref:
                    gpml_output += '>\n'
                    gpml_output += f'      <Xref identifier="{self.escape_xml(author.xref.identifier)}" dataSource="{self.escape_xml(author.xref.dataSource)}" />\n'
                    gpml_output += '    </Author>\n'
                else:
                    gpml_output += ' />\n'

            gpml_output += '  </Authors>\n'

        # Comments
        if pathway.comments:
            for comment in pathway.comments:
                if isinstance(comment, dict):
                    source = comment.get("source", "")
                    text = comment.get("text", comment.get("value", ""))
                else:
                    source = getattr(comment, 'source', '')
                    text = getattr(comment, 'value', '')

                if source:
                    gpml_output += f'  <Comment source="{self.escape_xml(source)}">{self.escape_xml(text)}</Comment>\n'
                else:
                    gpml_output += f'  <Comment>{self.escape_xml(text)}</Comment>\n'

        # Properties
        if pathway.properties:
            for prop in pathway.properties:
                if isinstance(prop, dict):
                    key = prop.get("key", "")
                    value = prop.get("value", "")
                else:
                    key = getattr(prop, 'key', '')
                    value = getattr(prop, 'value', '')

                gpml_output += f'  <Property key="{self.escape_xml(key)}" value="{self.escape_xml(value)}" />\n'

        # CitationRefs at pathway level
        if pathway.citationRefs:
            for citation_ref in pathway.citationRefs:
                ref_id = getattr(citation_ref, 'elementRef', '')
                gpml_output += f'  <CitationRef elementRef="{self.escape_xml(ref_id)}" />\n'

        # Graphics - board size
        if pathway.graphics:
            gpml_output += f'  <Graphics boardWidth="{pathway.graphics.boardWidth}" boardHeight="{pathway.graphics.boardHeight}"'
            if hasattr(pathway.graphics, 'backgroundColor') and pathway.graphics.backgroundColor:
                gpml_output += f' backgroundColor="{pathway.graphics.backgroundColor}"'
            gpml_output += ' />\n'

        # DataNodes
        if pathway.dataNodes:
            gpml_output += '  <DataNodes>\n'
            for datanode in pathway.dataNodes:
                gpml_output += self.write_datanode(datanode)
            gpml_output += '  </DataNodes>\n'


        # Interactions
        if pathway.interactions:
            gpml_output += '  <Interactions>\n'
            for interaction in pathway.interactions:
                gpml_output += self.write_interaction(interaction)
            gpml_output += '  </Interactions>\n'

        if hasattr(pathway, 'groups') and pathway.groups:
            gpml_output += '  <Groups>\n'
            for group in pathway.groups:
                gpml_output += self.write_group(group)
            gpml_output += '  </Groups>\n'

        # Citations
        if pathway.citations:
            gpml_output += '  <Citations>\n'
            for citation in pathway.citations:
                gpml_output += self.write_citation(citation)
            gpml_output += '  </Citations>\n'

        gpml_output += '</Pathway>\n'

        return gpml_output

    def write_citation(self, citation: Citation) -> str:
        """
        Write a Citation element to GPML format.

        Args:
            citation: Citation object

        Returns:
            str: GPML XML string for the Citation
        """
        element_id = self.escape_xml(citation.elementId)

        gpml_output = f'    <Citation elementId="{element_id}">\n'

        if citation.xref:
            identifier = self.escape_xml(citation.xref.identifier)
            data_source = self.escape_xml(citation.xref.dataSource)

            # Normalize data source names for GPML
            if data_source.lower() == "pubmed":
                data_source = "pubmed"
            elif data_source.lower() == "doi":
                data_source = "doi"

            gpml_output += f'      <Xref identifier="{identifier}" dataSource="{data_source}" />\n'

        if hasattr(citation, 'url') and citation.url:
            gpml_output += f'      <Url link="{self.escape_xml(citation.url.link)}" />\n'

        gpml_output += '    </Citation>\n'

        return gpml_output

        # Citations
        """"
        if pathway.citations:
            pass #skip for now
            gpml_output += '  <Citations>\n'
            for citation in pathway.citations:
                citation_id = citation.elementId
                # Remove 'citation_' prefix if present
                if citation_id.startswith('citation_'):
                    citation_id = citation_id[9:]

                gpml_output += f'    <Citation elementId="{self.escape_xml(citation_id)}">\n'

                if citation.xref:
                    gpml_output += f'      <Xref identifier="{self.escape_xml(citation.xref.identifier)}" dataSource="{self.escape_xml(citation.xref.dataSource)}" />\n'
                if hasattr(citation, 'url') and citation.url:
                    gpml_output += f'      <Url link="{self.escape_xml(citation.url.link)}" />\n'

                gpml_output += '    </Citation>\n'
            gpml_output += '  </Citations>\n'
        """

        gpml_output += '</Pathway>\n'

        return gpml_output