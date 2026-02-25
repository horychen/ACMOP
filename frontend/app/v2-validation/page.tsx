import StatorValidator from "@/components/StatorValidator";
import MotorWindingSimulator from "@/components/MotorWindingEngr";

export default function V2ValidationPage() {
    return (
        <div className="container mx-auto py-8 flex flex-col gap-8">
            <StatorValidator />
            <MotorWindingSimulator />
        </div>
    );
}
