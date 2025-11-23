import { DesignSpecs, PerformanceMetrics } from '../types';

// This is a mock service since we don't have the API key in this environment yet.
// In a real implementation, this would call the Google Generative AI API.

export const suggestSpecsFromDescription = async (description: string): Promise<Partial<DesignSpecs>> => {
    console.log('Mock AI Suggestion for:', description);
    // Return some dummy data based on keywords
    if (description.toLowerCase().includes('drone')) {
        return {
            ratedPower: 0.5,
            ratedSpeed: 5000,
            outerDiameterLimit: 80,
            axialLengthLimit: 60
        };
    }
    return {
        ratedPower: 5,
        ratedSpeed: 3000
    };
};

export const analyzeDesignResult = async (specs: DesignSpecs, performance: PerformanceMetrics): Promise<string> => {
    console.log('Mock AI Analysis for:', specs, performance);
    return `This design looks promising with an efficiency of ${performance.efficiency.toFixed(1)}%. The torque density is acceptable, but consider optimizing the slot depth for better thermal management.`;
};
