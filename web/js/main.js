document.getElementById('input').addEventListener('input', spaceMainInputChanged);
document.getElementById('input').addEventListener('input', transformMainInputChanged);
installInfoTooltipLayer();
updateTransformControls();
const initialWorkflow = location.hash.slice(1);
setWorkflow(WORKFLOW_CONFIG[initialWorkflow] ? initialWorkflow : 'polynomial', { keepInput: true, skipHash: true });
init();
