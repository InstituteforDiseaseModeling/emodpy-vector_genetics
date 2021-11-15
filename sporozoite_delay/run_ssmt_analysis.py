from analyzers.InsetChartAnalyzer import InsetChartAnalyzer
from analyzers.VectorGeneticsAnalyzer import VectorGeneticsAnalyzer
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis

if __name__ == "__main__":
    platform = Platform('CALCULON')
    analysis = PlatformAnalysis(platform=platform, experiment_ids=["4dac6575-df43-ec11-9ecd-9440c9bee941"],
                                analyzers=[InsetChartAnalyzer], # VectorGeneticsAnalyzer, InsetChartAnalyzer],
                                analyzers_args=[],
                                analysis_name="Sporozoite Delay",
                                extra_args=dict(partial_analyze_ok=True)
                                )

    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)