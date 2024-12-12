from analyzers.SummaryReportAnalyzer import SummaryReportAnalyzer
from analyzers.VectorStatsAnalyzerYearly import VectorStatsAnalyzerYearly
from analyzers.VectorGeneticsAnalyzer import VectorGeneticsAnalyzer
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis

if __name__ == "__main__":
    platform = Platform('CALCULON')
    analysis = PlatformAnalysis(platform=platform, experiment_ids=["8c22fca1-beb5-ef11-aa1a-b88303911bc1"],
                                analyzers=[SummaryReportAnalyzer, VectorGeneticsAnalyzer, VectorStatsAnalyzerYearly],
                                analyzers_args=[],
                                analysis_name="MIR drive high resistance high transmission",
                                extra_args=dict(partial_analyze_ok=True)
                                )

    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)