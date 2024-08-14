using  CSV, DataFrames


struct CobraProperties <: RobotProperties 
  minVol::Unitful.Volume
  maxVol::Unitful.Volume
  maxASP::Unitful.Volume
  positions::Vector{DeckPosition}
  compatible_stocks::Vector{DataType}
  CobraProperties(minVol,maxVol,maxASP,positions,compatible_stocks)=length(positions)==2 ? new(minVol,maxVol,maxASP,positions,compatible_stocks) : error("Cobra must have two defined deck positions")
end 

mutable struct CobraConfiguration <: RobotConfiguration
  source::String
  destination::String
  ASPRow::String
  ASPCol::Int
  washtime::Int # time in ms 
  ASPVol::Vector{Real}
  ASPPad::Real
  path::Vector{AbstractString}
  liquidclasses::Vector{AbstractString}
  pause::Bool
  predispenses::Int64
  cobra_path::String
end

struct Cobra <:Robot
  name::AbstractString 
  properties::CobraProperties 
  configuration::CobraConfiguration

end


struct SoftLinxSettings
  name::AbstractString
  n_loops::AbstractString
  cobrapath::AbstractString
end 

dwp96_2ml=JLIMS.Container("deep_well_plate_96_2_ml",2u"mL",(8,12))
dwp96_1ml=JLIMS.Container("deep_well_plate_96_1_ml",1u"mL",(8,12))
wp96=JLIMS.Container("plate_96",200u"µL",(8,12))
wp384=JLIMS.Container("plate_384",80u"µL",(16,24))

const cobra_names=Dict{JLIMS.Container,String}(
  dwp96_2ml=>"Deep Well 2 ml",
  dwp96_1ml=>"Deep Well - 1 ml",
  wp96=>"96 Costar",
  wp384=>"384 Well p/n 3575 3576")

const default_cobra_deck_1=DeckPosition("Deck 1",1,[
  dwp96_2ml,
  dwp96_1ml,
  wp384,
  wp96
])

const default_cobra_deck_2=DeckPosition("Deck 2",1,[
  dwp96_2ml,
  dwp96_1ml,
  wp384,
  wp96
])


cobra_default = Cobra("Default Cobra",
CobraProperties(0.3u"µL",40u"µL",750u"µL",[default_cobra_deck_1,default_cobra_deck_2],[JLIMS.LiquidStock]),
CobraConfiguration("N/A","N/A","A",1,5000,[0,0,0,0],1.1,["","","",""],["Water"],true,0,"C:\\Users\\Dell\\Dropbox (University of Michigan)\\JensenLab\\Cobra\\")
)

#############################################
# Wrap AcuteML Dependency into a single function 
#############################################


function fill_protocol_template(settings::CobraConfiguration)
  Source=settings.source
  Destination=settings.destination
  ASPRow=settings.ASPRow
  ASPCol=settings.ASPCol
  WashTime=settings.washtime
  ASPVol=settings.ASPVol
  Path=settings.path
  LiquidClass=settings.liquidclasses
  DispensePause=settings.usedispensepause
  PredispenseCount=settings.predispensecount


  template="""
<?xml version="1.0"?>
<Protocol xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <DispenseSpeed>49.6</DispenseSpeed>
  <DispenseOffset>0.25</DispenseOffset>
  <UseDispensePauseMode>$(DispensePause)</UseDispensePauseMode>
  <PredispenseCount>$(PredispenseCount)</PredispenseCount>
  <PredispenseLocation>Source</PredispenseLocation>
  <DeckLocations>
    <PlateSelection>
      <DeckNumber>1</DeckNumber>
      <AdapterName>None</AdapterName>
      <Plates>
        <PlateUseDescription>
          <Name>$(Source)</Name>
          <Use>Aspirate</Use>
        </PlateUseDescription>
      </Plates>
    </PlateSelection>
    <PlateSelection>
      <DeckNumber>2</DeckNumber>
      <AdapterName>None</AdapterName>
      <Plates>
        <PlateUseDescription>
          <Name>$(Destination)</Name>
          <Use>Dispense</Use>
        </PlateUseDescription>
      </Plates>
    </PlateSelection>
  </DeckLocations>
  <RepeatCount>1</RepeatCount>
  <ASPDistanceToBottom>2</ASPDistanceToBottom>
  <AspirateChannel1StartWell>
    <Row>$(ASPRow)</Row>
    <Column>$(ASPCol)</Column>
  </AspirateChannel1StartWell>
  <StagePostionForAspirate>1</StagePostionForAspirate>
  <AspirateIsSelected>true</AspirateIsSelected>
  <DispenseIsSelected>true</DispenseIsSelected>
  <DispensePattern>Serial</DispensePattern>
  <ParallelDispenseChannels>ChannelsAtoD</ParallelDispenseChannels>
  <WashIsSelected>true</WashIsSelected>
  <WashStrategy>ReservoirA</WashStrategy>
  <WashTimeAMsec>$(WashTime)</WashTimeAMsec>
  <WashTimeBMsec>0</WashTimeBMsec>
  <PurgeAfterDispense>false</PurgeAfterDispense>
  <DispenseChannels>
    <DispenseChannelInfo>
      <AspirateVolume>$(ASPVol[1])</AspirateVolume>
      <DispenseType>AirBacked</DispenseType>
      <DispenseVolume>0.1</DispenseVolume>
      <UseBottleCheck>false</UseBottleCheck>
      <UsesCustomFile>true</UsesCustomFile>
      <CustomFileName>$(Path[1])</CustomFileName>
      <LiquidClassName>$(LiquidClass[1])</LiquidClassName>
      <AspirateAdjustPercentage>0</AspirateAdjustPercentage>
      <DispenseAdjustPercentage>0</DispenseAdjustPercentage>
    </DispenseChannelInfo>
    <DispenseChannelInfo>
      <AspirateVolume>$(ASPVol[2])</AspirateVolume>
      <DispenseType>AirBacked</DispenseType>
      <DispenseVolume>0.2</DispenseVolume>
      <UseBottleCheck>false</UseBottleCheck>
      <UsesCustomFile>true</UsesCustomFile>
      <CustomFileName>$(Path[2])</CustomFileName>
      <LiquidClassName>$(LiquidClass[2])</LiquidClassName>
      <AspirateAdjustPercentage>0</AspirateAdjustPercentage>
      <DispenseAdjustPercentage>0</DispenseAdjustPercentage>
    </DispenseChannelInfo>
    <DispenseChannelInfo>
      <AspirateVolume>$(ASPVol[3])</AspirateVolume>
      <DispenseType>AirBacked</DispenseType>
      <DispenseVolume>0.3</DispenseVolume>
      <UseBottleCheck>false</UseBottleCheck>
      <UsesCustomFile>true</UsesCustomFile>
      <CustomFileName>$(Path[3])</CustomFileName>
      <LiquidClassName>$(LiquidClass[3])</LiquidClassName>
      <AspirateAdjustPercentage>0</AspirateAdjustPercentage>
      <DispenseAdjustPercentage>0</DispenseAdjustPercentage>
    </DispenseChannelInfo>
    <DispenseChannelInfo>
      <AspirateVolume>$(ASPVol[4])</AspirateVolume>
      <DispenseType>AirBacked</DispenseType>
      <DispenseVolume>0.5</DispenseVolume>
      <UseBottleCheck>false</UseBottleCheck>
      <UsesCustomFile>true</UsesCustomFile>
      <CustomFileName>$(Path[4])</CustomFileName>
      <LiquidClassName>$(LiquidClass[4])</LiquidClassName>
      <AspirateAdjustPercentage>0</AspirateAdjustPercentage>
      <DispenseAdjustPercentage>0</DispenseAdjustPercentage>
    </DispenseChannelInfo>
  </DispenseChannels>
  <PurgeAfterDispenseHeight>9</PurgeAfterDispenseHeight>
  <DispenseTestMode>None</DispenseTestMode>
</Protocol>
"""

    return template 

end 



function fill_softlinx_template(settings::SoftLinxSettings)


    name=settings.name
    n_loops=settings.n_loops
    filepath=settings.cobrapath

    template= """
    <Protocol mc:Ignorable="sap sap2010 sads" PreProtocolWizard="{x:Null}" ActivityLabel="" DisplayName="$(name)" HasConstraints="False" sap2010:WorkflowViewState.IdRef="Protocol_1" SLXId="94f2a419-8579-4634-bba3-e4cdb5050ddb" ToolTip="" UserComments="" isActive="True" isSetup="True"
 xmlns="clr-namespace:Hudson.Workflow.Activities;assembly=Hudson.Workflow.Activities"
 xmlns:hcc="clr-namespace:Hudson.Common.Communications;assembly=Hudson.Common"
 xmlns:hwab="clr-namespace:Hudson.Workflow.Activities.Base;assembly=SoftLinxBaseActivities"
 xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006"
 xmlns:p="http://schemas.microsoft.com/netfx/2009/xaml/activities"
 xmlns:s="clr-namespace:System;assembly=mscorlib"
 xmlns:sads="http://schemas.microsoft.com/netfx/2010/xaml/activities/debugger"
 xmlns:sap="http://schemas.microsoft.com/netfx/2009/xaml/activities/presentation"
 xmlns:sap2010="http://schemas.microsoft.com/netfx/2010/xaml/activities/presentation"
 xmlns:scg="clr-namespace:System.Collections.Generic;assembly=mscorlib"
 xmlns:x="http://schemas.microsoft.com/winfx/2006/xaml">
  <Protocol.Activities>
    <scg:List x:TypeArguments="p:Activity" Capacity="4">
      <LoopActivity ActivityLabel="" DisplayName="Loop Process" HasConstraints="False" sap2010:WorkflowViewState.IdRef="LoopActivity_1" SLXId="0859e942-9da6-494a-b0ac-7f2096172081" ToolTip="Loop $(n_loops) times." UserComments="loop through each protocol file in the folder" isActive="True" isSetup="True">
        <LoopActivity.Activities>
          <scg:List x:TypeArguments="p:Activity" Capacity="8">
            <ModifyVariableActivity Text="{x:Null}" CommandLine="dispense_file = path&amp;&quot;DispenseProtocol&quot;&amp;iLoop&amp;&quot;.xml&quot;" Description="" DisplayName="Modify Variable" HasConstraints="False" sap2010:WorkflowViewState.IdRef="ModifyVariableActivity_2" SLXId="4dd62ca7-9f97-4540-af3a-42a0567673d6" ToolTip="" UserComments="rename the dispense file  according to its iLoop number " isActive="True" isCanceled="False" isSetup="True">
              <ModifyVariableActivity.Arguments>
                <ModifyVariableActivityArguments DisplayExpression="path&amp;&quot;DispenseProtocol&quot;&amp;iLoop&amp;&quot;.xml&quot;" DisplayVariable="dispense_file" Expression="path&amp;&quot;DispenseProtocol&quot;&amp;iLoop&amp;&quot;.xml&quot;" VariableName="dispense_file" />
              </ModifyVariableActivity.Arguments>
              <ModifyVariableActivity.TimeConstraints>
                <scg:List x:TypeArguments="hwab:TimeConstraint" Capacity="0" />
              </ModifyVariableActivity.TimeConstraints>
            </ModifyVariableActivity>
            <AdvancedInstrumentActivity CommandLine="Run" Description="Arguments: dispense_file" DisplayName="Advanced Cobra" HasConstraints="False" sap2010:WorkflowViewState.IdRef="AdvancedInstrumentActivity_1" SLXId="431dcce0-f40a-4460-99c3-ac6e02468696" ToolTip="Command: Run&#xA;Arguments: dispense_file" UserComments="" isActive="True" isCanceled="False" isSetup="True">
              <AdvancedInstrumentActivity.Arguments>
                <InstrumentActivityArguments Address="{x:Reference __ReferenceID8}" ResultVariable="{x:Null}" AddinType="Cobra" Command="Run">
                  <InstrumentActivityArguments.Arguments>
                    <scg:List x:TypeArguments="x:Object" Capacity="4">
                      <x:String>dispense_file</x:String>
                    </scg:List>
                  </InstrumentActivityArguments.Arguments>
                  <InstrumentActivityArguments.hWnd>
                    <s:IntPtr />
                  </InstrumentActivityArguments.hWnd>
                </InstrumentActivityArguments>
              </AdvancedInstrumentActivity.Arguments>
              <AdvancedInstrumentActivity.TimeConstraints>
                <scg:List x:TypeArguments="hwab:TimeConstraint" Capacity="0" />
              </AdvancedInstrumentActivity.TimeConstraints>
            </AdvancedInstrumentActivity>
            <ModifyVariableActivity Text="{x:Null}" CommandLine="iLoop = iLoop +1" Description="" DisplayName="Modify Variable" HasConstraints="False" sap2010:WorkflowViewState.IdRef="ModifyVariableActivity_1" SLXId="b2a9e9ea-a7b6-4327-bd6e-8fc370ab7dac" ToolTip="" UserComments="iterate loop" isActive="True" isCanceled="False" isSetup="True">
              <ModifyVariableActivity.Arguments>
                <ModifyVariableActivityArguments DisplayExpression="iLoop +1" DisplayVariable="iLoop" Expression="iLoop +1" VariableName="iLoop" />
              </ModifyVariableActivity.Arguments>
              <ModifyVariableActivity.TimeConstraints>
                <scg:List x:TypeArguments="hwab:TimeConstraint" Capacity="0" />
              </ModifyVariableActivity.TimeConstraints>
            </ModifyVariableActivity>
          </scg:List>
        </LoopActivity.Activities>
        <LoopActivity.Activities2>
          <scg:List x:TypeArguments="p:Activity" Capacity="0" />
        </LoopActivity.Activities2>
        <LoopActivity.Arguments>
          <LoopActivityArguments CountExpression="$(n_loops)" IsBoolean="False" IsCounted="True" IsInfinite="False" TrueExpression="" />
        </LoopActivity.Arguments>
        <LoopActivity.TimeConstraints>
          <scg:List x:TypeArguments="hwab:TimeConstraint" Capacity="0" />
        </LoopActivity.TimeConstraints>
      </LoopActivity>
    </scg:List>
  </Protocol.Activities>
  <Protocol.Activities2>
    <scg:List x:TypeArguments="p:Activity" Capacity="0" />
  </Protocol.Activities2>
  <Protocol.InitialValues>
    <scg:Dictionary x:TypeArguments="x:String, hwab:Variable" />
  </Protocol.InitialValues>
  <Protocol.Interfaces>
    <hwab:Interface x:Key="{x:Reference __ReferenceID4}" SetupData="{x:Array Type=x:String}" AddinType="Cobra">
      <hwab:Interface.Address>
        <hcc:SLAddress x:Name="__ReferenceID4" Name="Cobra" Workcell="SoftLinx" />
      </hwab:Interface.Address>
    </hwab:Interface>
    <hwab:Interface x:Key="{x:Reference __ReferenceID5}" SetupData="{x:Array Type=x:String}" AddinType="Files">
      <hwab:Interface.Address>
        <hcc:SLAddress x:Name="__ReferenceID5" Name="Files" Workcell="SoftLinx" />
      </hwab:Interface.Address>
    </hwab:Interface>
    <hwab:Interface x:Key="{x:Reference __ReferenceID6}" SetupData="{x:Null}" AddinType="FileTemplateEditor">
      <hwab:Interface.Address>
        <hcc:SLAddress x:Name="__ReferenceID6" Name="FileTemplateEditor" Workcell="SoftLinx" />
      </hwab:Interface.Address>
    </hwab:Interface>
    <hwab:Interface x:Key="{x:Reference __ReferenceID7}" SetupData="{x:Null}" AddinType="CSVReader">
      <hwab:Interface.Address>
        <hcc:SLAddress x:Name="__ReferenceID7" Name="CSVReader" Workcell="SoftLinx" />
      </hwab:Interface.Address>
    </hwab:Interface>
  </Protocol.Interfaces>
  <Protocol.TimeConstraints>
    <scg:List x:TypeArguments="hwab:TimeConstraint" Capacity="0" />
  </Protocol.TimeConstraints>
  <Protocol.Variables>
    <hwab:VariableList SLXHost="{x:Null}">
      <hwab:Variable x:TypeArguments="hwab:Interface" Value="{x:Reference __ReferenceID0}" x:Key="SoftLinx.Files" Name="SoftLinx.Files" Prompt="False">
        <hwab:Variable.Default>
          <hwab:Interface SetupData="{x:Null}" x:Name="__ReferenceID0" AddinType="Files">
            <hwab:Interface.Address>
              <hcc:SLAddress Name="Files" Workcell="SoftLinx" />
            </hwab:Interface.Address>
          </hwab:Interface>
        </hwab:Variable.Default>
      </hwab:Variable>
      <hwab:Variable x:TypeArguments="hwab:Interface" Value="{x:Reference __ReferenceID1}" x:Key="SoftLinx.FileTemplateEditor" Name="SoftLinx.FileTemplateEditor" Prompt="False">
        <hwab:Variable.Default>
          <hwab:Interface SetupData="{x:Null}" x:Name="__ReferenceID1" AddinType="FileTemplateEditor">
            <hwab:Interface.Address>
              <hcc:SLAddress Name="FileTemplateEditor" Workcell="SoftLinx" />
            </hwab:Interface.Address>
          </hwab:Interface>
        </hwab:Variable.Default>
      </hwab:Variable>
      <hwab:Variable x:TypeArguments="hwab:Interface" Value="{x:Reference __ReferenceID2}" x:Key="SoftLinx.CSVReader" Name="SoftLinx.CSVReader" Prompt="False">
        <hwab:Variable.Default>
          <hwab:Interface SetupData="{x:Null}" x:Name="__ReferenceID2" AddinType="CSVReader">
            <hwab:Interface.Address>
              <hcc:SLAddress Name="CSVReader" Workcell="SoftLinx" />
            </hwab:Interface.Address>
          </hwab:Interface>
        </hwab:Variable.Default>
      </hwab:Variable>
      <hwab:Variable x:TypeArguments="x:String" x:Key="dispense_file" Default="&quot;&quot;" Name="dispense_file" Prompt="False" Value="&quot;&quot;" />
      <hwab:Variable x:TypeArguments="x:Double" x:Key="iLoop" Default="1" Name="iLoop" Prompt="False" Value="1" />
      <hwab:Variable x:TypeArguments="hwab:Interface" Value="{x:Reference __ReferenceID3}" x:Key="SoftLinx.Cobra" Name="SoftLinx.Cobra" Prompt="False">
        <hwab:Variable.Default>
          <hwab:Interface SetupData="{x:Array Type=x:String}" x:Name="__ReferenceID3" AddinType="Cobra">
            <hwab:Interface.Address>
              <hcc:SLAddress x:Name="__ReferenceID8" Name="Cobra" Workcell="SoftLinx" />
            </hwab:Interface.Address>
          </hwab:Interface>
        </hwab:Variable.Default>
      </hwab:Variable>
      <hwab:Variable x:TypeArguments="x:String" x:Key="path" Default="$(filepath)" Name="path" Prompt="False" Value="$(filepath)" />
    </hwab:VariableList>
  </Protocol.Variables>
  <sap2010:WorkflowViewState.ViewStateManager>
    <sap2010:ViewStateManager>
      <sap2010:ViewStateData Id="ModifyVariableActivity_2" sap:VirtualizedContainerService.HintSize="413,79" />
      <sap2010:ViewStateData Id="AdvancedInstrumentActivity_1" sap:VirtualizedContainerService.HintSize="413,79" />
      <sap2010:ViewStateData Id="ModifyVariableActivity_1" sap:VirtualizedContainerService.HintSize="413,79" />
      <sap2010:ViewStateData Id="LoopActivity_1" sap:VirtualizedContainerService.HintSize="453,461" />
      <sap2010:ViewStateData Id="Protocol_1" sap:VirtualizedContainerService.HintSize="453,597" />
    </sap2010:ViewStateManager>
  </sap2010:WorkflowViewState.ViewStateManager>
  <sads:DebugSymbol.Symbol>dztDOlxVc2Vyc1xEZWxsXERvY3VtZW50c1wyMDIzXzA4XzA4X0plbnNlbkRpc3BlbnNlQ29icmEuc2x2cAUBAZUBDAEBDwc+FgEFEg0ZJgEOGg0qKgELKw0yJgEI</sads:DebugSymbol.Symbol>
</Protocol>

    """

return template

end 



########## Cobra Dispense Files ###########

function cobraCSV(design::DataFrame,source::JLIMS.Container,destination::JLIMS.Container,path::String)
    
  r_in,c_in=source.shape
  n_in=r_in*c_in
  r_out,c_out=destination.shape
  n_out=r_out*c_out
  if mod(ncol(design),4) != 0
      error("number of reagent columns must be a multiple of 4. Pad the design with empty reagent columns if needed")
  end 

  if ncol(design) > n_in 
      error("number of design reagents exceeds the size of the source plate. Split the design across multiple source plates")
  end 

  if nrow(design) != n_out
      error("output plate size does not match design. Expected design to have $(n_out) rows. Pad design with empty experiments if needed or split the design across multiple destination plates.")
  end 

  if !isdir(path)
      mkdir(path)
  end 

  filenames=String[]
  for j in 1:ncol(design)

      x=design[:,j]
      x=reshape(x,r_out,c_out)
      x=DataFrame(x,:auto)
      filename=joinpath(path,"DispenseFile$(j).csv")
      filenames[j]=filename
      CSV.write(filename,x,;writeheader=false)
  end 
  return filenames
end 


function fill_design(design::DataFrame,source::JLIMS.Container,destination::JLIMS.Container)
  x=Matrix(design)

  
  r_in,c_in=source.shape
  n_in=r_in*c_in
  r_out,c_out=destination.shape
  n_out=r_out*c_out

  y=zeros(n_out,n_in)

  a,b=size(x)

  y[1:a,1:b].=x
  return DataFrame(y,:auto)
end 


function snake_order(r,c;channel=isodd) # generate the order of dispenses for a plate of any size. The cobra snakes through each plate: ie. starts in row 1 and goes across the row forward, shifts down to row 2 and goes backward across the row  etc. 
  n=r*c
  N=1:n
  N=reshape(N,r,c)
  idxs=[0 for _ in 1:n]
  counter=1
  for i in 1:r
    for j in (channel(i) ? (1:c) : reverse(1:c))
      idxs[counter]=N[i,j]
      counter+=1
    end 
  end 
  return idxs 
end


function design2protocols(directory::String,design::DataFrame,source::JLIMS.Container,destination::JLIMS.Container,robot::Cobra)
  pad=robot.configuration.ASPPad
  maxASP=robot.properties.maxASP
  maxShot=robot.properties.maxVol
  cobra_location=robot.configuration.cobra_path
  # Check for issues with the design

    r_in,c_in=source.shape
    n_in=r_in*c_in
    r_out,c_out=destination.shape
    n_out=r_out*c_out
    des=Matrix(design)
    r,c=size(des)
    if r != n_in
      error("output plate size does not match design. Expected design to have $(n_out) rows. Pad design with empty experiments if needed or split the design across multiple destination plates.")
    end 
    if r != length(config.liquidclasses)
      error("each of the $c reagents must have a specified liquid class")
    end 
    if c > n_out 
      error("number of design reagents exceeds the size of the source plate. Split the design across multiple source plates")
    end 
    if mod(r,4) != 0 
        ArgumentError("The number of design columns must be a multiple of 4")
    end 
    filecounter=1
    rows=string.(collect('A':4:'Z'))[1:fld(r_in,4)]
    cols=1:c_in
    pos=collect(Iterators.product(rows,cols))
    positions=reshape(pos,prod(size(pos)),1)
    protocols=String[]
    configs=Cobra[]
  # Start the protocol generation process  
    for set in 1:fld(c,4) # loop through sets of four locations in the source plate 
      idxs=4*(set-1)+1:4*(set-1)+4
      d = Matrix(des[:,idxs]) # grab the appropriate source wells from the design 
      while sum(d)>0 # while there is still liquid to be dispensed in the set 
        to_dispense=zeros(size(d))
        fileidxs=Int64[]
        for channel in 1:4
          disp_order=snake_order(r_out,c_out; channel= (isodd(channel) ? iseven : isodd)) # the even channels snake opposite to the odd ones ( yes i know this is confusing, I suggest you watch the cobra run)
          total_vols=sum(to_dispense[:,channel])*pad
          for idx in disp_order
            if  total_vols < maxASP  # grab all of the rows we can complete with a single pass 
                margin= max((maxASP - total_vols),0)
                to_dispense[idx,channel]=min(d[idx,channel],maxShot,margin/pad)
                total_vols=sum(to_dispense[:,channel])*pad
            else
              break 
            end 
          end 
          x=reshape(to_dispense[:,channel],r_out,c_out)
          x=DataFrame(x,:auto)
          filename=joinpath(directory,"DispenseFile$(filecounter).csv")
          push!(fileidxs,filecounter)
          filecounter+=1
          CSV.write(filename,x)   
        end 
        d=d.-to_dispense

        cobra_name=robot.name
        cobra_properties=robot.properties
        
        cobra_config=deepcopy(robot.configuration)

        cobra_config.ASPRow=positions[set][1]
        cobra_config.ASPCol=positions[set][2]
        cobra_config.ASPVol=vec(sum(to_dispense,dims=1))*pad
        cobra_config=config.liquidclasses[idxs]
        cobra_config.source=cobra_names[source]
        cobra_config.destination=cobra_names[destination]

        cobra_config.path=[cobra_location*"$(experiment_name)\\DispenseFile$(k).csv" for k in fileidxs]
        cobra_config.cobra_path=cobra_location
        configured_cobra=Cobra(cobra_name,cobra_properties,cobra_config)
        push!(configs,configured_cobra)
        #settings=CobraSettings(source,destination,positions[set][1],positions[set][2],washtime,vec(sum(to_dispense,dims=1))*pad,cobra_path,liquidclasses[idxs],dispensepause,predispensecount)
        push!(protocols,fill_protocol_template(settings))
      end  
    end 
    return protocols,configs
end 
        


function protocols2softlinx(n_protocols::Int64,experiment_name::String,cobra_location::String)
    settings=SoftLinxSettings(experiment_name,string(n_protocols),string(cobra_location,experiment_name,"\\"))
    softlinx_out=fill_softlinx_template(settings)
    return softlinx_out
end 

# plate to plate mixer
function cobra_mixer(sources::Vector{JLIMS.Stock},destinations::Vector{JLIMS.Stock},robot::Cobra;quiet=false,timelimit=10)
  source_deck=robot.properties.positions[1]
  destination_deck=robot.properties.positions[2]
  any(map(x -> !in(x.well.container,source_deck.compatible_containers)||!in(typeof(x),robot.properties.compatible_stocks),sources)) ? nothing : error("Robot $(robot.name) is incapable of using at least one of the requested source stocks")
  if any(map(x->!in(typeof(x),robot.properties.compatible_stocks)||!in(x.well.container,destination_deck.compatible_containers),destinations)) 
      error("Robot $(robot.name) is incapable of making at least one of the requested stocks")
  end 


  source_labware=unique(map(x->x.well.labwareid,sources))

  sl_idx=map(x->findall(y->y.well.labwareid==x,sources),source_labware)
  source_containers=map(x->sources[x[1]].well.container,sl_idx)

  destination_labware=unique(map(x->x.well.labwareid,desitnations))
  dl_idx=map(x->findall(y->y.well.labwareid==x,destinations),destination_labware)
  destination_containers=map(x->destiantions[x[1]].well.container,dl_idx)
  sl=length(source_labware)
  dl=length(destination_labware)


  concentrations=stock_concentration_array(sources)
  target_quantities=stock_quantity_array(destinations)
  concentrations=concentrations[:,DataFrames.names(target_quantities)]
      
  concentrations=Float64.(ustrip.(Matrix(concentrations)))

  target_quantities=Float64.(ustrip.(Matrix(target_quantities)))

  source_quantities= map(x->x.quantity,available_sources) 
  source_units= preferred_stock_quantity(available_sources)
  source_quants=uconvert.(source_quantities,source_units)
      
  S,I=size(concentrations)
  D,I=size(target_quantities)
  model=Model(Gurobi.Optimizer)
  if quiet 
      set_silent(model)
  end 
      set_attribute(model,"TimeLimit",timelimit)

      
      @variable(model,q[1:S,1:D] >=0)
      @variable(model,lw[1:sl,1:dl],Bin)
      @variable(model,multipass_well[1:S,1:D],Bin)
      @variable(model, multipass_stock[1:S],Bin)
      @constraint(model, q'*concentrations .== target_quantities) # find the quantity of each stock needed to make the target 

      for s in 1:S 
        @constraint(model, sum(q[s,:])*robot.configuration.ASPPad <= stock_quants[s])
      end 

      for s in 1:sl 
        for d in 1:dl 
          @constraint(model, !lw[s,d]=>{sum(q[sl_idx[s],dl_idx[d]]) == 0})
        end 
      end 
      for s in 1:S
        for d in 1:D
          @constraint(model, !multipass_well[s,d]=> {q[s,d] <=robot.properties.maxVol}) # flag wells that need multiple passes from a source
        end 
      end 
      for s in 1:S
        @constraint(model, !multipass_stock[s]=> {sum(q[s,:]) <= robot.properties.maxASP}) # flag sources that need multiple passes 
      end 

      @objective(model, Min,sum(lw)) # minimize the number of labware pairings needed (number of protocols to be written) 
      optimize!(model)

      w =JuMP.Value.(lw)
      @constraint(model, lw .==w)
      set_objective_sense(model, MOI.FEASIBILITY_SENSE)
      @objective(model, Min,sum(mulitpass_well)) # minimize the number of wells that need multiple passes to fill (dispense volume greater than max shot volume of cobra)
      optimize!(model)

      m =JuMP.Value.(mutlipass_well)
      @constraint(model, multipass_well .==m)
      set_objective_sense(model, MOI.FEASIBILITY_SENSE)
      @objective(model, Min,sum(multipass_stock)) # minimize the number of revisits to a particular stock (sum(dispense volumes) greater than max aspiration volume for cobra)
      optimize!(model)

      quants=JuMP.value.(q)
  

      
      quantities=Any[]
      for stock in available_sources
          un=preferred_stock_quantity(stock)
          if stock in srcs 
              idx=findfirst(x->x==stock,srcs)
              val=quants[idx]*un
              if isa(val,Unitful.Mass)
                  if 0u"g" <val < robot.minMass
                      error("The $val mass required for the stock in well #$(stock.well.id) is too small for a $(typeof(robot)) to accurately transfer")
                  elseif  val > robot.maxMass
                      error("The $val mass required for the stock in well #$(stock.well.id) is too large for a $(typeof(robot)) to accurately transfer")
                  end 
              elseif isa(val,Unitful.Volume)
                  if 0u"L" <val < robot.minVol
                      error("The $val volume required for the stock in well #$(stock.well.id) is too small for a $(typeof(robot)) to accurately transfer")
                  elseif  val > robot.maxVol
                      error("The $val volume required for the stock in well #$(stock.well.id) is too large for a $(typeof(robot)) to accurately transfer")
                  end 
              end 
              push!(quantities,quants[idx]*un) # commit the transfer to the transfers vector for this destination stock 
              
          else
              push!(quantities,0*un) # commit a 0 to the transfers vector
          end 
      end 
      wellnames=["Well$(i)" for i in map(x->x.well.id,destinations)]
      transfers=DataFrame(quantities,wellnames) # volumes in µL
      tf=ustrip.(transfers) 
      n_pairs=sum(w) 
      protocols=Tuple{DataFrame,Vector{JLIMS.Stock},Vector{JLIMS.Stock}}[]
      for i in 1:sl 
        for j in 1:dl 
          if w[i,j]==1 
            if i ==j 
              error("cannot schedule the same labware in two positions at the same time")
            end 
            active_sources=sources[sl_idx[i]]
            active_destinations=destinations[dl_idx[j]]
            rows=prod(source_containers[i].shape)
            cols=prod(destination_containers[j].shape)
            df = DataFrame(zeros(rows,cols))
            row_idxs=map(x->x.well.wellindex,active_sources)
            col_idxs=map(x->x.well.wellindex,active_destinations)
            df[row_idxs,col_idxs] .= tf[sl_idx[i],dl_idx[j]]
            push!(protocols,(df,active_sources,active_destinations))
          end 
        end 
      end 

      transfer_table=DataFrame(Source=Integer[],Destination=Integer[],Quantity=Real[],Unit=AbstractString[])
      r=nrow(transfers)
      c=ncol(transfers)
      for col in 1:c
          for row in 1:r 
              val=transfers[row,col]
              quantity=ustrip(val)
              if quantity==0 
                  continue 
              else 
                  source=sources[row].well.id 
                  destination=destinations[col].well.id
                  un=string(unit(val))
                  push!(transfer_table,(source,destination,quantity,un))
              end 
          end 
      end 

      
      return transfer_table, protocols 
  end 
 




"""
    cobra(design::DataFrame,directory::String,source::String,destination::String, liquidclasses::Vector{String} ;pad::Real=1.1,pause::Bool=true,predispensecount::Int=0,washtime::Int=5000, kwargs...)

Create Cobra dipsense instructions for microplate source to destination operations

  ## Arguments 
  * `design`: a (# of destinations) x (# of sources) dataframe containing the volume of each dispense in µL.
  * `directory`: the ouput directory of the files. directory will automatically be created if it doesn't already exist  
  * `source`: The source plate type. See `keys(containers)` for available options.
  * `destination`: The destination plate type. See `keys(containers)` for available options. 
  * `liquidclasses`: A vector of liquid class types. There must be a class for each reagents. Use "Water" as default. 

  ## Keyword Arguments 
  * `pad=1.1`: Pad the aspriation volume by a multiplier to ensure that the stock doesn't run out during dispense. Default is 1.1 (10%).
  * `maxASP=800`: The maximum aspiration volume for a channel. We have 1 ml syringes, but we only allow them to aspirate up to 800 µl. 
  * `maxShot=40`: The maximum shot volume for a dispense in µl. 
  * `pause=true`: When true, Cobra pauses over each well during dispense; when false, it moves continuously.
  * `washtime=8000`: The length of the wash step in milliseconds between each dispensing operation. Default is 8000 ms.
  * `predispensecount=0`: Number of pre-dispense shots Cobra performs before dispensing 
  *  `cobra_location="C:\\Users\\Dell\\Dropbox (University of Michigan)\\JensenLab\\Cobra\\"`: The location of the cobra's working directory on the driving machine.
"""
function cobra_instructor(directory::String,protocol_name::AbstractString,design::DataFrame,sources::Vector{JLIMS.Stock},destinations::Vector{JLIMS.Stock},robot::Cobra;liquidclass="Water")
    allequal(map(x->x.well.labwareid,sources)) ? nothing : error("All source stocks must be on the same labware")
    allequal(map(x.well.labwareid,destinations)) ? nothing : error("All destination stocks must be on the same labware")
    source=sources[1].well.container
    destination=destinations[1].well.container 
    robot.configuration.liquidclasses=[liquidclass for _ in 1:ncol(design)] #hard code the same liquid class for each reagent 
    if ~isdir(directory)
      mkdir(directory)
    end 

    full_dir=joinpath(directory,protocol_name)
    if ~isdir(full_dir)
      mkdir(full_dir)
    end 
    protocols,cobra_configs=design2protocols(full_dir,design,source,destination,robot)
    protocol_names=["DispenseProtocol$(k).xml" for k in 1:length(protocols)]
    full_protocol_names=joinpath.((full_dir,),protocol_names)
    map((loc,protocol) -> write(loc,protocol),full_protocol_names,protocols)
    n_protocols=length(protocols)
    softlinx=protocols2softlinx(n_protocols,protocol_name,robot.configuration.cobra_path)
    full_softlinx_name=joinpath(full_dir,"$(protocol_name).slvp")
    write(full_softlinx_name,softlinx)
    write(joinpath(full_dir,"config.json"),JSON.json(cobra_configs))
    #print("entire $(experiment_name) folder must be moved to Dropbox -> JensenLab -> Cobra")
end




function cobra(directory::AbstractString,sources::Vector{JLIMS.Stock},destinations::Vector{JLIMS.Stock},robot::Cobra;kwargs...)
  tt, protocols=cobra_mixer(sources,destinations,robot;kwargs...)
  names=random_protocol_name(length(protocols))
  for i in eachindex(protocols)
    cobra_instructor(directory,names[i],protocols[i][1],protocols[i][2],protocols[i][3],robot;kwargs...)
  end 
  return tt 
end 








#=
settings=cobra_settings()
dispensefiles=["./src/Cobra/2023_06_13_test_cobra/test1/Dispensefile$(k).csv" for k in 1:96]
liquidclasses=["Water" for _ in 1:96]
templatefile= "./src/Cobra/JensenCobraTemplate.xml"


protocols=dispensefiles2protocols(dispensefiles,templatefileliquidclasses)

testxmls = ["./src/Cobra/xml_test/test$(i).xml" for i in 1:length(protocols)]
map((loc,protocol) -> write(loc,protocol),testxmls,protocols)

print(String(read("./src/Cobra/JensenDispenseSoftlinxTemplate.xml")))

=#

#=
test_design=CSV.read("./src/Cobra/2023_06_13_test_cobra/test1.csv",DataFrame)

source="96-1mL-deepwell"
destination="384-well"
filepath="./src/Cobra/CobraDispenseTest/"
liquidclasses=["Water" for _ in 1:96]
CobraDispense(test_design,filepath,source,destination,liquidclasses)
=#

#= 
using JLDispense,DataFrames
design=ones(384,96)
design=DataFrame(design,:auto)
lqs=["Water" for _ in 1:96]
cobra(design, "/Users/BDavid/Desktop/cobra_test/",containers[:dWP96_2ml],containers[:WP384],lqs)
=#