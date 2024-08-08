#include "MeshKernelApi/MKStateUndoAction.hpp"

std::unique_ptr<meshkernelapi::MKStateUndoAction> meshkernelapi::MKStateUndoAction::Create(MeshKernelState& mkState)
{
    return std::make_unique<MKStateUndoAction>(mkState);
}

meshkernelapi::MKStateUndoAction::MKStateUndoAction(MeshKernelState& mkState) : m_mkStateReference(mkState)
{
    // Copy the pointers
    m_mkState.m_mesh1d = mkState.m_mesh1d;
    m_mkState.m_network1d = mkState.m_network1d;
    m_mkState.m_mesh2d = mkState.m_mesh2d;
    m_mkState.m_contacts = mkState.m_contacts;
    m_mkState.m_curvilinearGrid = mkState.m_curvilinearGrid;
    m_mkState.m_meshOrthogonalization = mkState.m_meshOrthogonalization;
    m_mkState.m_curvilinearGridFromSplines = mkState.m_curvilinearGridFromSplines;
    m_mkState.m_curvilinearGridOrthogonalization = mkState.m_curvilinearGridOrthogonalization;
    m_mkState.m_curvilinearGridLineShift = mkState.m_curvilinearGridLineShift;
    m_mkState.m_projection = mkState.m_projection;
    m_mkState.m_state = mkState.m_state;
}

void meshkernelapi::MKStateUndoAction::SwapContents()
{
    std::swap(m_mkState.m_mesh1d, m_mkStateReference.m_mesh1d);
    std::swap(m_mkState.m_mesh2d, m_mkStateReference.m_mesh2d);
    std::swap(m_mkState.m_network1d, m_mkStateReference.m_network1d);
    std::swap(m_mkState.m_contacts, m_mkStateReference.m_contacts);
    std::swap(m_mkState.m_curvilinearGrid, m_mkStateReference.m_curvilinearGrid);
    std::swap(m_mkState.m_state, m_mkStateReference.m_state);
}

void meshkernelapi::MKStateUndoAction::DoCommit()
{
    SwapContents();
}

void meshkernelapi::MKStateUndoAction::DoRestore()
{
    SwapContents();
}
