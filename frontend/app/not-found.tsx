export default function NotFound() {
    return (
        <div className="flex items-center justify-center h-screen">
            <div className="text-center">
                <h2 className="text-2xl font-bold mb-2">404 - Page Not Found</h2>
                <p className="text-muted-foreground">The page you're looking for doesn't exist.</p>
            </div>
        </div>
    )
}

export const dynamic = 'force-dynamic'
