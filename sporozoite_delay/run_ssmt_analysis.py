from analyzers.InsetChartAnalyzer import InsetChartAnalyzer
from analyzers.VectorGeneticsAnalyzer import VectorGeneticsAnalyzer
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis

if __name__ == "__main__":
    platform = Platform('CALCULON')
    analysis = PlatformAnalysis(platform=platform, experiment_ids=["8bbe4c48-6326-ec11-9ecd-9440c9bee941"],
                                analyzers=[VectorGeneticsAnalyzer], # VectorGeneticsAnalyzer, InsetChartAnalyzer],
                                analyzers_args=[],
                                analysis_name="Sporozoite Delay VG",
                                extra_args=dict(partial_analyze_ok=True)
                                )

    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)