/*=========================================================================

  Program:   Visualization Library
  Module:    DS2SPtsF.cc
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

Description:
---------------------------------------------------------------------------
This file is part of the Visualization Library. No part of this file
or its contents may be copied, reproduced or altered in any way
without the express written consent of the authors.

Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen 1993, 1994 

=========================================================================*/
#include "DS2SPtsF.hh"

void vlDataSetToStructuredPointsFilter::Update()
{
  vlDataSetFilter::Update();
}

void vlDataSetToStructuredPointsFilter::PrintSelf(ostream& os, vlIndent indent)
{
  if (this->ShouldIPrint(vlDataSetToStructuredPointsFilter::GetClassName()))
    {
    this->PrintWatchOn(); // watch for multiple inheritance
    
    vlStructuredPoints::PrintSelf(os,indent);
    vlDataSetFilter::PrintSelf(os,indent);
    
    this->PrintWatchOff(); // stop worrying about it now
    }
}
