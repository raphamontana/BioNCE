import { Step, Stepper, StepLabel } from '@material-ui/core';
import DatasetComponent from "../components/ml/DatasetComponent"
import ModelSelectComponent from "../components/ml/ModelSelectComponent"
import ResultsComponent from "../components/ml/ResultsComponent"
import ReviewComponent from "../components/ml/ReviewComponent"
import Layout from '../components/layout/Layout';

const ML = () => {
  const [dataset, setDataset] = React.useState([]);
  const [model, setModel] = React.useState('');
  const [activeStep, setActiveStep] = React.useState(0);
  const steps = ['Set dataset', 'Select model', 'Run prediction', "Results"];

  const handleNext = () => {
    setActiveStep((prevActiveStep) => prevActiveStep + 1);
  };

  const handleBack = () => {
    setActiveStep((prevActiveStep) => prevActiveStep - 1);
  };

  const handleReset = () => {
    setDataset([]);
    setModel('');
    setActiveStep(0);
  };

  return (
    <Layout>
      <Stepper activeStep={activeStep}>
        {steps.map((label, index) => {
          const stepProps = {};
          const labelProps = {};
          return(
            <Step key={label} {...stepProps}>
              <StepLabel {...labelProps}>{label}</StepLabel>
            </Step>
          );
        })}
      </Stepper>
      {activeStep === 0 ? (
        <DatasetComponent dataset={dataset} setDataset={setDataset} handleNext={handleNext} />
      ) : activeStep === 1 ? (
        <ModelSelectComponent model={model} setModel={setModel} handleBack={handleBack} handleNext={handleNext} />
      ) : activeStep === 2 ? (
        <ReviewComponent dataset={dataset} model={model} handleBack={handleBack} handleNext={handleNext} />
      ) : (
        <ResultsComponent handleReset={handleReset} />
      )}
    </Layout>
  );
};

export default ML;