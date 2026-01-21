from sc_pipeline import SCPipeline

datasets = [
    {
        'dataset_name': 'Adam',
        'input_path': '数据集/Adam/data.h5', 
        'output_base': 'Results',
        'target_col': 'cell_ontology_class',
        'root_cell_type': 'ureteral cell',
        'organism': 'Mouse',
        
        # 实验模式开关
        'use_experiment': False, 
        # 'experiment_files': {'embeddings': 'path/to/emb.txt', 'pred': 'path/to/pred.txt'},
        
        # 可选：手动指定基因 (如无则自动)
        # 'specific_genes': {
        #     'Ureteric bud': ['Rps26', 'Ptma', 'Lsp1', 'Rbp1'],
        #     'PT': ['Lrp2', 'Aldob', 'Slc34a1', 'Gatm', 'Kap'],
        #     'CM': ['Col1a2', 'Col3a1', 'Tmsb4x', 'Gucy1a3'],
        #     'Distal': ['Slc12a3', 'Calb1', 'Wnk1', 'Emcn'],
        #     'Endothelial': ['Kdr', 'Flt1', 'Pecam1', 'Ets1'],
        #     'Podocytes': ['Nphs1', 'Nphs2', 'Wt1', 'Synpo', 'Podxl'],
        #     'Loop of Henle': ['Umod', 'Slc12a1', 'Spp1'],
        #     'Stromal': ['Dcn', 'Mgp', 'Sparc']
        # }
    }
]

def main():
    for config in datasets:
        print(f"\n{'='*40}")
        print(f"🎬 开始处理: {config['dataset_name']}")
        print(f"{'='*40}")
        
        # inspect_data 参数设置
        # inspect_data=False 关闭数据检查，否则开启数据检查
        # inspect_data=Ture 打开数据检查，查看h5数据集是否规范，是否包含了测量矩阵，基因信息等关键数据
        pipe = SCPipeline(config, inspect_data=False)
        pipe.load_data()
        pipe.run_preprocessing()
        pipe.find_markers()
        pipe.plot_markers() 
        
        try:
            pipe.run_enrichment() # 这里会自动调用 plot_kegg_network
        except Exception as e:
            print(f"⚠️ 富集分析跳过: {e}")
            
        pipe.run_trajectory()
        pipe.save_data()
        
    print("\n🎉 完成!")

if __name__ == "__main__":
    main()