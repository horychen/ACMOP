import StatorValidator from "@/components/StatorValidator";
import MotorWindingSimulator from "@/components/MotorWindingEngr";
import AnimatedSvgVisualizer from "@/components/AnimatedSvgVisualizer";

export default function V2ValidationPage() {
    return (
        <div className="container mx-auto py-8 flex flex-col gap-8">
            <h1 className="text-3xl font-bold mb-2">V2 Design Validation</h1>
            <AnimatedSvgVisualizer />
            <StatorValidator />
            <MotorWindingSimulator />
        </div>
    );
}
