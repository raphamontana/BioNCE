import { Step, Stepper, StepLabel } from '@material-ui/core';
import DatasetComponent from "../components/ml/DatasetComponent"
import ModelSelectComponent from "../components/ml/ModelSelectComponent"
import ResultsComponent from "../components/ml/ResultsComponent"
import ReviewComponent from "../components/ml/ReviewComponent"
import Layout from '../components/layout/Layout';

const ML = () => {
  const [activeStep, setActiveStep] = React.useState(0);
  const [skipped, setSkipped] = React.useState(new Set());
  const steps = ['Set dataset', 'Select model', 'Run prediction', "Results"];

  const isStepSkipped = (step) => {
    return skipped.has(step);
  };

  const handleNext = () => {
    let newSkipped = skipped;
    if (isStepSkipped(activeStep)) {
      newSkipped = new Set(newSkipped.values());
      newSkipped.delete(activeStep);
    }
    setActiveStep((prevActiveStep) => prevActiveStep + 1);
    setSkipped(newSkipped);
  };

  const handleBack = () => {
    setActiveStep((prevActiveStep) => prevActiveStep - 1);
  };

  const handleReset = () => {
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
        <DatasetComponent handleNext={handleNext} />
      ) : activeStep === 1 ? (
        <ModelSelectComponent handleBack={handleBack} handleNext={handleNext} />
      ) : activeStep === 2 ? (
        <ReviewComponent handleBack={handleBack} handleNext={handleNext} />
      ) : (
        <ResultsComponent handleReset={handleReset} />
      )}
    </Layout>
  );
};

export default ML;