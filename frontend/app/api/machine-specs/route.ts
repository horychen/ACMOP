import { NextResponse } from 'next/server';
import fs from 'fs';
import path from 'path';

export async function GET() {
    try {
        // Adjust this path to match your local environment structure
        // Assuming frontend is at c:\_Codes\ACMOP\frontend and backend is at c:\_Codes\ACMOP\backend
        const filePath = path.join(process.cwd(), '../backend/codes4/machine_specifications.json');

        if (!fs.existsSync(filePath)) {
            return NextResponse.json({ error: 'Specifications file not found' }, { status: 404 });
        }

        const fileContent = fs.readFileSync(filePath, 'utf-8');
        const data = JSON.parse(fileContent);

        return NextResponse.json(data);
    } catch (error) {
        console.error('Error reading machine specs:', error);
        return NextResponse.json({ error: 'Failed to load specifications' }, { status: 500 });
    }
}
