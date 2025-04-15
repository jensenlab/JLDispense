using  CSV, DataFrames


struct SoftLinxSettings
  name::AbstractString
  n_loops::AbstractString
  cobrapath::AbstractString
end 

struct CobraProtocolSettings
  source::AbstractString
  destination::AbstractString
  ASPRow::AbstractString
  ASPCol::Int
  ASPVol::Vector{Real}
  path::Vector{AbstractString}
  liquidclasses::Vector{AbstractString}
end



# define cobra configuration 

struct CobraSettings <: InstrumentSettings 
  pause::Bool
  predispenses::Int64
  cobra_path::AbstractString # path to dispense files on the machine running the cobra 
  washtime::Int # time in ms 
  ASPPad::Real # a multiplicataive padding factor
  ASPDistance::Real # distance in mm
end




struct CobraHead <: TransferHead
  nozzles::AbstractArray{Nozzle}
end 

struct CobraDeckPosition <: DeckPosition 
  labware::Set{Type{<:Labware}}
end 


CobraConfiguration = Configuration{CobraHead,Deck{CobraDeckPosition},CobraSettings} 

const cobra_nozzle = ContinuousNozzle(0.3u"µL",40u"µL",800u"µL",20u"µL",1.1,true)


const cobra_head = CobraHead(fill(cobra_nozzle,4,1))

const cobra_compat_labware= Set([JLConstants.DeepWP96,JLConstants.WP96,JLConstants.WP384])

const cobra_position_1=CobraDeckPosition(cobra_compat_labware)
const cobra_position_2=CobraDeckPosition(cobra_compat_labware)

const cobra_deck = [cobra_position_1,cobra_position_2]

const cobra_settings= CobraSettings(true,0,"C:\\Users\\Dell\\Dropbox (University of Michigan)\\JensenLab\\Cobra\\",8000,1.1,2)

const cobra =CobraConfiguration(cobra_head,cobra_deck,cobra_settings)


liquidclass(::JLIMS.Stock) = "Water" 
 

### define deck access functions 


function can_aspirate(h::CobraHead, d::CobraDeckPosition,l::Labware) 
  return can_place(l,d)
end
function can_dispense(h::CobraHead,d::CobraDeckPosition,l::Labware) 
  return can_place(l,d)
end

# define masks 

function masks(h::CobraHead,l::JLConstants.WellPlate) # for generic 96 well plates, we will define a separate method for 384 well plates 
  C= length(nozzles(h))
  Wi,Wj=shape(l)
  Pi,Pj = 5,12
  W = falses(Wi,Wj)
  P=falses(Pi,Pj)
  function Ma(w::Integer,p::Integer,c::Integer) 
      # w=wells, p=positions, c=channels
      1 <= w <= Wi*Wj || return false 
      1 <= p <= Pi*Pj || return false 
      1 <= c <= C || return false 
      wi,wj=cartesian(W,w)
      pm,pn=cartesian(P,p)
      return wi == c+pm-1 && wj == pn 
  end 
  C= length(nozzles(h))
  Wi,Wj=shape(l)
  Pi,Pj = 11,12
  W = falses(Wi,Wj)
  P=falses(Pi,Pj)
  function Md(w::Integer,p::Integer,c::Integer) 
      # w=wells, p=positions, c=channels
      1 <= w <= Wi*Wj || return false 
      1 <= p <= Pi*Pj || return false 
      1 <= c <= C || return false 
      wi,wj=cartesian(W,w)
      pm,pn=cartesian(P,p)
      return wi == c+pm-4 && wj == pn 
  end 
  return Ma,Md
end 

function masks(h::CobraHead,l::JLConstants.WP384) 
  C= length(nozzles(h))
  Wi,Wj=shape(l)
  Pi,Pj = 10,12
  W = falses(Wi,Wj)
  P=falses(Pi,Pj)
  function Ma(w::Integer,p::Integer,c::Integer) 
      # w=wells, p=positions, c=channels
      1 <= w <= Wi*Wj || return false 
      1 <= p <= Pi*Pj || return false 
      1 <= c <= C || return false 
      wi,wj=cartesian(W,w)
      pm,pn=cartesian(P,p)
      return wi == 2*(c-1)+pm && wj == pn 
  end 
  C= length(nozzles(h))
  Wi,Wj=shape(l)
  Pi,Pj = 22,12
  W = falses(Wi,Wj)
  P=falses(Pi,Pj)
  function Md(w::Integer,p::Integer,c::Integer) 
      # w=wells, p=positions, c=channels
      1 <= w <= Wi*Wj || return false 
      1 <= p <= Pi*Pj || return false 
      1 <= c <= C || return false 
      wi,wj=cartesian(W,w)
      pm,pn=cartesian(P,p)
      return wi == 2*(c-1)+pm-6 && wj == pn 
  end 
  return Ma,Md
end 






const cobra_names=Dict{Type{<:Labware},AbstractString}(
  JLConstants.DeepWP96=>"Deep Well 2 ml",
  WP96=>"96 Costar",
  WP384=>"384 Well p/n 3575 3576")













#############################################
# Wrap AcuteML Dependency into a single function 
#############################################


function fill_protocol_template(settings::CobraSettings,protocol::CobraProtocolSettings)
  Source=protocol.source
  Destination=protocol.destination
  ASPRow=protocol.ASPRow
  ASPCol=protocol.ASPCol
  WashTime=settings.washtime
  ASPVol=protocol.ASPVol
  Path=protocol.path
  LiquidClass=protocol.liquidclasses
  DispensePause=settings.pause
  PredispenseCount=settings.predispenses
  distance=settings.ASPDistance


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
  <ASPDistanceToBottom>$(distance)</ASPDistanceToBottom>
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
#=
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

=#
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


function design2protocols(directory::AbstractString,design::DataFrame,source::JLIMS.Labware,destination::JLIMS.Labware,config::CobraConfiguration)
  pad=settings(config).ASPPad
  maxASP=ustrip.(uconvert.(u"µL",map(x->x.maxASP,nozzles(head(config))) ))
  maxShot=ustrip(uconvert(u"µL",map(x->x.maxVol,nozzles(head(config)))))
  cobra_location=settings(config).cobra_path
  # Check for issues with the design

    r_in,c_in=shape(source)
    n_in=r_in*c_in
    r_out,c_out=shape(destination)
    n_out=r_out*c_out
    des=Matrix(design)
    r,c=size(des)
    if r != n_in
      error("output plate size does not match design. Expected design to have $(n_out) rows. Pad design with empty experiments if needed or split the design across multiple destination plates.")
    end 
    if r != length(robot.configuration.liquidclasses)
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
            if  total_vols < maxASP[channel]  # grab all of the rows we can complete with a single pass 
                margin= max((maxASP[channel] - total_vols),0)
                to_dispense[idx,channel]=min(d[idx,channel],maxShot[channel],margin/pad)
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
          CSV.write(filename,x;writeheader=false)   
        end 
        d=d.-to_dispense



        ASPRow=positions[set][1]
        ASPCol=positions[set][2]
        ASPVol=vec(sum(to_dispense,dims=1))*pad
        liquidclasses=fill("Water",4)
        source=cobra_names[source]
        destination=cobra_names[destination]
        protocol_name=basename(directory)
        path=[cobra_location*"$(protocol_name)\\DispenseFile$(k).csv" for k in fileidxs]
        protocol_settings=CobraProtocolSettings(source,destination,ASPRow,ASPCol,ASPVol,path,liquidclasses)
        push!(configs,config)
        #settings=CobraSettings(source,destination,positions[set][1],positions[set][2],washtime,vec(sum(to_dispense,dims=1))*pad,cobra_path,liquidclasses[idxs],dispensepause,predispensecount)
        push!(protocols,fill_protocol_template(protocol_settings,config))
      end  
    end 
    return protocols,configs
end 
        


function protocols2softlinx(n_protocols::Int64,experiment_name::String,cobra_location::String)
    settings=SoftLinxSettings(experiment_name,string(n_protocols),string(cobra_location,experiment_name,"\\"))
    softlinx_out=fill_softlinx_template(settings)
    return softlinx_out
end 





"""
    dispense(config::CobraConfiguration, design::DataFrame, directory::AbstractString,protocol_name::AbstractString,source::Labware,destination::Labware)

Create Cobra dipsense instructions for well plate dispensing 

  ## Arguments 
  * `config`: A Configuration object defining the cobra 
  * `design`: a (# of destinations) x (# of sources) dataframe containing the volume of each dispense in µL.
  * `directory`: the ouput directory of the files. directory will automatically be created if it doesn't already exist  
  * `source`: The source plate type. See `keys(containers)` for available options.
  * `destination`: The destination plate type. See `keys(containers)` for available options. 
  * `liquidclasses`: A vector of liquid class types. There must be a class for each reagents. Use "Water" as default. 


"""
function dispense(config::CobraConfiguration, design::DataFrame, directory::AbstractString,protocol_name::AbstractString,source::Labware,destination::Labware) 
    full_dir=joinpath(directory,protocol_name)
    if ~isdir(full_dir)
      mkdir(full_dir)
    end 
    protocols,cobra_configs=design2protocols(full_dir,dispenses,source,destination,config)
    protocol_names=["DispenseProtocol$(k).xml" for k in 1:length(protocols)]
    full_protocol_names=joinpath.((full_dir,),protocol_names)
    map((loc,protocol) -> write(loc,protocol),full_protocol_names,protocols)
    n_protocols=length(protocols)
    softlinx=protocols2softlinx(n_protocols,protocol_name,settings(config).cobra_path)
    full_softlinx_name=joinpath(full_dir,"$(protocol_name).slvp")
    write(full_softlinx_name,softlinx)
    write(joinpath(full_dir,"config.json"),JSON.json(cobra_configs))
    #print("entire $(experiment_name) folder must be moved to Dropbox -> JensenLab -> Cobra")
    
    return nothing 
    # include the wasted sources because we pad the aspiration  Well # 1 is designated for waste. 
end

function dispense(config::CobraConfiguration,design::DataFrame,directory::AbstractString,protocol_name::AbstractString,source::Type{<:Labware},destination::Type{<:Labware})
  src=source(1,"source_instance")
  dst=destinatino(1,"destination_instance")
  return dispense(config,design,directory,protocol_name,src,dst)
end 



#=
function mixer(directory::AbstractString,sources::Vector{T},destinations::Vector{U},robot::Cobra;kwargs...) where  {T <: JLIMS.Stock,U <:JLIMS.Stock}
  design=dispense_solver(sources,destinations,robot,minimize_overdrafts!,minimize_transfers!;kwargs...)

  srcs_needed=filter(x->ustrip(sum(design[x,:])) > 0,eachindex(sources))

  srcs=sources[srcs_needed]

  des=design[srcs_needed,:]
  dest_labware=unique(map(x->x.well.labwareid,destinations))
  dest_idxs=[findall(x->x.well.labwareid==i,destinations) for i in dest_labware]
  src_labware=unique(map(x->x.well.labwareid,srcs))
  src_idxs=[findall(x->x.well.labwareid==i,srcs) for i in src_labware]
  tt=DataFrame[]
  pn=AbstractString[]
  for i in eachindex(src_labware)
    for j in eachindex(dest_labware)
        dispenses=des[src_idxs[i],dest_idxs[j]]
        if sum(Matrix(dispenses)) > 0u"µL"
            protocol_name=random_protocol_name()
            transfer_table=cobra(directory,protocol_name,dispenses,srcs[src_idxs[i]],destinations[dest_idxs[j]],robot;kwargs...)
            push!(pn,protocol_name)
            push!(tt,transfer_table)
        end 
    end 
  end 
  return pn,tt
end 
=#







#=
using JLDispense,JLD2

sources=JLD2.load("./src/Mixer/dwp_stocks.jld2")["stocks"]
destinations=JLD2.load("./src/Mixer/stock_targets.jld2")["out"]
alts=JLD2.load("./src/Mixer/example_stocks.jld2")["stocks"]

protocol_name,transfer_table=mixer("/Users/BDavid/Desktop/",sources,destinations,cobra_default)

=# 