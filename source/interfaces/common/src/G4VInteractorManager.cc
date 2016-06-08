//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VInteractorManager.cc,v 1.11 2001/12/03 08:07:45 barrand Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// G.Barrand

#include <stdlib.h>
#include <string.h>

#include "g4std/algorithm"

#include "G4VInteractorManager.hh"

#define NewString(str)  \
 ((str) != NULL ? (strcpy((char*)malloc((unsigned)strlen(str) + 1), str)) : NULL)

/***************************************************************************/
G4VInteractorManager::G4VInteractorManager (
)
:argc(0)
,argv(NULL)
,mainInteractor(NULL)
,secondaryLoopEnabled(TRUE)
,alreadyInSecondaryLoop(FALSE)
,exitSecondaryLoop(0)
,parentInteractor(NULL)
,createdInteractor(NULL)
,creationString(NULL)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
}
/***************************************************************************/
G4VInteractorManager::~G4VInteractorManager (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(argv!=NULL) {
    for(G4int argi=0;argi<argc;argi++) {
      if(argv[argi]!=NULL) free(argv[argi]);
    }
    free (argv);
  }
  argv = NULL;
  argc = 0;
  dispatchers.clear();
  preActions.clear();
  postActions.clear();
  shells.clear();
  secondaryLoopEnabled = TRUE;
  alreadyInSecondaryLoop = FALSE;
  exitSecondaryLoop = 0;
}
/***************************************************************************/
void G4VInteractorManager::SetArguments (
 G4int  a_argc
,char** a_argv
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  // Free previous values.
  if(argv!=NULL) {
    for(G4int argi=0;argi<argc;argi++) {
      if(argv[argi]!=NULL) free(argv[argi]);
    }
    free(argv);
  }
  argv = NULL;
  argc = 0;
  // Set new values.
  if(a_argc!=0) {
    argv = (char**)malloc(a_argc * sizeof(char*));
    if(argv!=NULL) {
      argc = a_argc;
      for(G4int argi=0;argi<a_argc;argi++) {
	argv[argi] = (char*)NewString (a_argv[argi]);
      }
    }
  }
}
/***************************************************************************/
char** G4VInteractorManager::GetArguments (
 G4int* a_argc
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_argc!=NULL) *a_argc = argc;
  return argv;
}
/***************************************************************************/
void G4VInteractorManager::SetMainInteractor (
 G4Interactor a_main
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  mainInteractor = a_main;
}
/***************************************************************************/
G4Interactor G4VInteractorManager::GetMainInteractor (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return mainInteractor;
}
/***************************************************************************/
void G4VInteractorManager::EnableSecondaryLoop (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  secondaryLoopEnabled = TRUE;
}
/***************************************************************************/
void G4VInteractorManager::DisableSecondaryLoop (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  secondaryLoopEnabled = FALSE;
}
/***************************************************************************/
void G4VInteractorManager::AddDispatcher (
 G4DispatchFunction a_dispatcher
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_dispatcher==NULL) return;
  if(G4std::find(dispatchers.begin(),dispatchers.end(),a_dispatcher)!=dispatchers.end()) return;
  dispatchers.push_back(a_dispatcher);
}
/***************************************************************************/
void G4VInteractorManager::RemoveDispatcher (
 G4DispatchFunction a_dispatcher
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4std::vector<G4DispatchFunction>::iterator it;
  for (it = dispatchers.begin(); it != dispatchers.end(); it++) {
    if (*it == a_dispatcher) {
      dispatchers.erase(it);
      break;
    }
  }
}
/***************************************************************************/
void G4VInteractorManager::DispatchEvent (
 void* a_event
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4int dispatchern = dispatchers.size();
  G4DispatchFunction func;
  for(G4int count=0;count<dispatchern;count++) {
    func = dispatchers[count];
    if(func!=NULL) {
      if(func(a_event)==true) return;
    }
  }
}
/***************************************************************************/
void G4VInteractorManager::AddSecondaryLoopPreAction (
 G4SecondaryLoopAction a_preAction
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_preAction==NULL) return;
  if(G4std::find(preActions.begin(),preActions.end(),a_preAction)!=preActions.end()) return;
  preActions.push_back(a_preAction);
}
/***************************************************************************/
void G4VInteractorManager::SecondaryLoopPreActions (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4int preActionn = preActions.size();
  for(G4int count=0;count<preActionn;count++) {
    if(preActions[count]!=NULL) preActions[count]();
  }
}
/***************************************************************************/
void G4VInteractorManager::AddSecondaryLoopPostAction (
 G4SecondaryLoopAction a_postAction
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_postAction==NULL) return;
  if(G4std::find(postActions.begin(),postActions.end(),a_postAction)!=postActions.end()) return;
  postActions.push_back(a_postAction);
}
/***************************************************************************/
void G4VInteractorManager::SecondaryLoopPostActions (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  G4int postActionn = postActions.size();
  for(G4int count=0;count<postActionn;count++) {
    if(postActions[count]!=NULL) postActions[count]();
  }
}
/***************************************************************************/
void G4VInteractorManager::SecondaryLoop (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(Inited()==FALSE) return;

  if(secondaryLoopEnabled==FALSE) return;
  
  if (alreadyInSecondaryLoop==FALSE) {
    G4cout << "------------------------------------------" << G4endl;
    G4cout << "You have entered a viewer secondary X event loop." << G4endl;
    G4cout << "Quit it with an 'Escape' viewer button" << G4endl;
    alreadyInSecondaryLoop   = TRUE;
    exitSecondaryLoop        = 0;
    SecondaryLoopPreActions  ();
    //for(G4int count=0;count<shelln;count++) XWidgetUniconify(shells[count]);
    void*                    event;
    while(1) {
      event = GetEvent();
      if(event==NULL) break;
      DispatchEvent  (event);
      if(exitSecondaryLoop!=0) break;
    }
    G4cout << "Secondary X event loop exited." << G4endl;
    SecondaryLoopPostActions ();
    }
}
/***************************************************************************/
void G4VInteractorManager::RequireExitSecondaryLoop (
 G4int a_code
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(secondaryLoopEnabled==FALSE) return;
  if(a_code==0)            a_code = 1;
  exitSecondaryLoop        = a_code;
  alreadyInSecondaryLoop   = FALSE;
  // for(G4int count=0;count<shelln;count++) XWidgetIconify(shells[count]);
  // if(shelln!=0)            XSync(XtDisplay(topWidget),False);
}
/***************************************************************************/
G4int G4VInteractorManager::GetExitSecondaryLoopCode (
) 
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return exitSecondaryLoop;
}
/***************************************************************************/
void G4VInteractorManager::AddShell (
 G4Interactor a_shell
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  if(a_shell==NULL) return;
  if(G4std::find(shells.begin(),shells.end(),a_shell)!=shells.end()) return;
  shells.push_back(a_shell);
}
/***************************************************************************/
void G4VInteractorManager::RemoveShell (
 G4Interactor a_shell
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{  
  G4std::vector<G4Interactor>::iterator it;
  for (it = shells.begin(); it != shells.end(); it++) {
    if (*it == a_shell) {
      shells.erase(it);
      break;
    }
  }
}
/***************************************************************************/
void G4VInteractorManager::SetParentInteractor (
 G4Interactor a_interactor
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  parentInteractor = a_interactor;
}
/***************************************************************************/
G4Interactor G4VInteractorManager::GetParentInteractor (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return parentInteractor;
}
/***************************************************************************/
void G4VInteractorManager::SetCreatedInteractor (
 G4Interactor a_interactor
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  createdInteractor = a_interactor;
}
/***************************************************************************/
G4Interactor G4VInteractorManager::GetCreatedInteractor (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return createdInteractor;
}
/***************************************************************************/
void G4VInteractorManager::SetCreationString (
 char* a_string
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  creationString = a_string;
}
/***************************************************************************/
char* G4VInteractorManager::GetCreationString (
)
/***************************************************************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
  return creationString;
}
