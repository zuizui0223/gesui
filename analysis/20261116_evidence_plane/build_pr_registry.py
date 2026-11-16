#!/usr/bin/env python3
"""Exhaustive real GitHub PR research registry on a frozen 2026-11-16 snapshot.

Every PR has a title, source branch, head commit SHA, state and file-path list.
Only already reviewed PRs receive scientific verdicts; unknown stays blank.
Main/staged inclusion is defined by identical Git blob SHA, never by README
mentions or coincidental filenames. Positive + null/STOP evidence preserved.
"""
from __future__ import annotations
import argparse,csv,json,os,time
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import requests

API="https://api.github.com"
REPO="zuizui0223/gesui"
SNAPSHOT=Path("provenance/PR_ALL_CURRENT_METADATA_SOURCE_SNAPSHOT_20261116.json")
HEADERS={"Accept":"application/vnd.github+json",
         "X-GitHub-Api-Version":"2022-11-28",
         "User-Agent":"gesui-source-pr-registry/20261116"}
# Deliberately blank for all other PRs; no result guessed from PR titles.
REVIEW={
14:("肯定＋制限","EEA生物観測局×2022","上流WWTP年別報告負荷の底生EQR地域外予測"),
16:("肯定・否定混合／欠測","EEA生物観測局×2022","欠測・農業・源負荷の対抗検査"),
17:("肯定＋地点依存","EEA生物観測局×2018/2022","2018年の再現と別地点・年別転送"),
81:("肯定＋制限","US NRSA UID×調査期/HUC8","実測物理環境の3期EPT外挿"),
82:("肯定＋感度","US NRSA UID×調査期/HUC8","200個体希薄化・分類群数"),
83:("肯定＋感度","US NRSA UID×調査期/HUC8","EPT FFG/HABIT/PTV"),
90:("追加予測不確実","EU H7×期間","報告源N/Pの追加H6予測なし"),
92:("方向対照で不確実","EU H7×期間","上流/下流施設負荷方向比較"),
93:("方向対照で不確実","EU H7×期間","下流調整後上流強度の予測不確実"),
95:("肯定＋得点感度","EU H7×期間","蝶Lycaenidae条件付き正"),
96:("層間残差は不確実","EU H7×期間","水EQR・植物・蝶の残差同調なし"),
97:("間接経路不確実","US NRSA同一訪問/HUC8","TN→Chla正、間接積は不確実"),
98:("肯定＋感度","US NRSA UID×3調査期/HUC8","河床埋没度とEPT固着形質"),
101:("肯定＋科別注意","EU H7×期間","蝶HNS科内種ラベル入替"),
102:("肯定＋時間null","EU H7×期間","植物NUTRI・蝶HNS地理正"),
103:("肯定＋共通6形質STOP","EU H7×期間","NUTRI対照5形質を個別調整"),
104:("肯定＋水質未確認","EU H7内異なるH7/真H6","植物Nと蝶HNSのH6内対応"),
105:("属間正・属内不確実","EU H7内異なるH7/真H6","植物Nの属間属内分解"),
106:("main統合済み","研究原本単位別","98元Git blobとnull/STOP保存"),
107:("水質媒介source STOP","EEA生物観測局×年/H7","同局9 H7期間/近傍1km 16上限"),
108:("方向特異性は不確実","EU H7×期間","WWTP上流強度の下流対照後上乗せなし"),
109:("方向維持＋外挿感度","EU H7×期間","固定記録数・都市化との対抗検査"),
110:("空間正＋時間外挿不確実","EU H7×期間","EQR-植物N逆関連、H6外不確実")
}
def api(url,token):
    head=dict(HEADERS)
    if token:head["Authorization"]="Bearer "+token
    for k in range(4):
        try:
            r=requests.get(url,headers=head,timeout=(8,25))
            if r.status_code==200:return r.json()
            if r.status_code in (403,429,500,502,503,504):
                time.sleep(2**k*2);continue
            raise RuntimeError("Official GitHub API "+str(r.status_code)+" "+url)
        except requests.RequestException:
            if k==3:raise
            time.sleep(2**k)
    raise RuntimeError("Official GitHub API retry exhausted "+url)
def change_files(n,token):
    out=[]
    for page in range(1,15):
        data=api(f"{API}/repos/{REPO}/pulls/{n}/files?per_page=100&page={page}",token)
        if not isinstance(data,list):raise RuntimeError("Non-list PR changed files "+str(n))
        out.extend(data)
        if len(data)<100:return out
    raise RuntimeError("Unexpected huge changed-file PR "+str(n))
def tree(ref,token):
    obj=api(f"{API}/repos/{REPO}/git/trees/{ref}?recursive=1",token)
    if obj.get("truncated") or not isinstance(obj.get("tree"),list):
        raise RuntimeError("Git tree incomplete "+ref)
    return {r["path"]:r["sha"] for r in obj["tree"] if r.get("type")=="blob"}
def is_result_doc(path):
    return path.startswith("docs/RESULTS_") or (path.startswith("results/") and path.rsplit("/",1)[-1] in ("RESULTS.md","RESULTS.txt","REPORT.md"))

def main(args):
    original=json.loads(SNAPSHOT.read_text())
    prs=original["prs"]
    if original.get("schema")!="GESUI_20261116_CURRENT_GITHUB_ALL_PR_METADATA_RAW_SNAPSHOT_V1" or len(prs)<100:
        raise RuntimeError("Complete real GitHub PR snapshot absent")
    if len({p["number"] for p in prs})!=len(prs):raise RuntimeError("Duplicate PR ID in source")
    if not args.token:raise RuntimeError("Only authorized public GitHub API read token needed")
    mainblobs=tree("main",args.token)
    staged=tree(args.branch,args.token)
    files={}
    with ThreadPoolExecutor(max_workers=8) as exe:
        fs={exe.submit(change_files,p["number"],args.token):p["number"] for p in prs}
        for future in as_completed(fs):
            files[fs[future]]=future.result()
    rows=[];doc_manifest=[];all_docs=set();notinmain=set();all_staged=set()
    for item in sorted(prs,key=lambda q:q["number"]):
        n=item["number"]
        changed=files[n]
        docs=[x for x in changed if x.get("status")!="removed" and is_result_doc(x["filename"])]
        paths=sorted({x["filename"] for x in docs})
        if len(paths)!=len(docs):raise RuntimeError("Repeated result doc path within PR "+str(n))
        inside=[x["filename"] for x in docs if mainblobs.get(x["filename"])==x.get("sha")]
        staging=[x["filename"] for x in docs if staged.get(x["filename"])==x.get("sha")]
        missing=[x["filename"] for x in docs if mainblobs.get(x["filename"])!=x.get("sha")]
        all_docs.update(paths);notinmain.update(missing);all_staged.update(staging)
        verdict,grain,conclusion=REVIEW.get(n,("","",""))
        rows.append({
          "pr_number":n,"PR_URL":item["html_url"],"title":item["title"],
          "head_branch":item["head_ref"],"head_commit_SHA":item["head_sha"],
          "base_branch":item["base_ref"],"state":item["state"],
          "is_merged":bool(item.get("merged_at")),"is_draft":bool(item["draft"]),
          "created_at":item["created_at"],"updated_at":item["updated_at"],
          "analysis_grain_reviewed":grain,"verdict_reviewed_only":verdict,
          "summary_reviewed_only":conclusion,
          "result_document_paths":"; ".join(paths),
          "result_document_count":len(paths),
          "source_identical_result_documents_in_main":"; ".join(inside),
          "source_identical_result_documents_staged":"; ".join(staging),
          "result_documents_missing_or_changed_in_main":"; ".join(missing),
          "identical_docs_main_count":len(inside),
          "identical_docs_staged_count":len(staging),
        })
        doc_manifest.append({
          "PR":n,"original_head_sha_snapshot":item["head_sha"],
          "original_changed_files_n":len(changed),
          "changed_original_result_docs":[{"path":x["filename"],"source_blob_SHA1":x.get("sha")} for x in docs],
          "verdict_unreviewed_intentionally_blank":n not in REVIEW,
        })
    args.output.mkdir(parents=True,exist_ok=True)
    dest=args.output/"PR_ALL_RESEARCH_REGISTRY_20261116.csv"
    with dest.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    manifest={"schema":"GESUI_REAL_ALL_PR_RESULTS_DOC_MANIFEST_V1","snapshot":str(SNAPSHOT),
      "all_PRs":doc_manifest,
      "current_PR_changed_file_list_retrieved_after_head_SHA_snapshot":True,
      "snapshot_head_SHA_may_have_changed_after_snapshot":True}
    (args.output/"PR_SOURCE_FILE_DOC_MANIFEST_20261116.json").write_text(
        json.dumps(manifest,ensure_ascii=False,indent=2)+"\n")
    summary={"schema":"GESUI_REAL_ALL_PR_RESULT_DOC_REGISTRY_AND_MAIN_GAPS_V1",
        "num_PRs_with_actual_original_metadata":len(prs),
        "num_PR_changed_file_lists_read":len(files),
        "known_verdicts":sum(bool(x["verdict_reviewed_only"]) for x in rows),
        "unreviewed_verdicts_left_blank":sum(not bool(x["verdict_reviewed_only"]) for x in rows),
        "different_PR_result_document_paths":len(all_docs),
        "different_PR_result_document_paths_not_exact_in_main":len(notinmain),
        "different_PR_result_document_paths_exact_in_staged_branch":len(all_staged),
        "GitHub_main_file_sha_exact_matching_only":True,
        "no_title_based_scientific_verdict_inference":True,
        "not_all_historical_PR_source_files_moved_into_main":True,
        "original_PR_head_snapshot_SHA_preserved":True}
    (args.output/"PR_REGISTRY_SUMMARY_20261116.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    doc=[
      "# 全PRの検証結果登録簿（2026-11-16）","",
      "本登録簿はGithubの全PRの元head SHAと実際の変更ファイル一覧を読み、単位・判定・元結果文書・mainに同一Git blobが存在するかを分けて記録したもの。",
      "",
      "判定が空欄なのは未精査であり、肯定・否定・停止のいずれとも推測していない。main READMEに言及があっても結果ファイルの統合とは数えない。",
      "",
      "- PR全数: **"+str(len(prs))+"**。",
      "- 結論精査済み: **"+str(summary["known_verdicts"])+"**、未判定: **"+str(summary["unreviewed_verdicts_left_blank"])+"**。",
      "- PR結果文書の異なるパス: **"+str(len(all_docs))+"**、mainに元blobが一致しない結果文書パス: **"+str(len(notinmain))+"**。",
      "",
      "全PRの行別テーブル: [原本PR登録簿CSV](../provenance/PR_ALL_RESEARCH_REGISTRY_20261116.csv)。詳細な元PR変更文書とGit blob SHA: [doc manifest JSON](../provenance/PR_SOURCE_FILE_DOC_MANIFEST_20261116.json)。",
      "",
      "## 主な精査済みPR",
      "",
      "| PR | 単位 | 判定 | head SHA（短縮） | 元結果文書のmain同一数 |",
      "|---|---|---|---|---|"
    ]
    for row in rows:
        if not row["verdict_reviewed_only"]:continue
        doc.append("| [#"+str(row["pr_number"])+"]("+row["PR_URL"]+") | "+row["analysis_grain_reviewed"]+
          " | "+row["verdict_reviewed_only"]+" | "+row["head_commit_SHA"][:12]+" | "+
          str(row["identical_docs_main_count"])+"/"+str(row["result_document_count"])+" |")
    doc.extend(["","## 同じ平面で比較する際の単位","",
      "EUは真正HydroBASINS H7×元の2期間の454行（有効植物NUTRI451行）。水生EQR・植物NUTRI・蝶HNS・WWTP報告N/Pを同じ流域と期間へ揃える。EEA実測負荷と底生EQRは元観測局からH7に再集計するが独立観測数は増えない。",
      "米国NRSAの3期はUID×survey-cycle、地域外検証はHUC8単位。欧州のモデルWWTP地理Wと米国の河床埋没度XEMBEDは別の曝露であり同一の生態学的効果量としてプールしない。",
      "",
      "GitHubにPRのhead SHAを記録しても、元ブランチ削除後の恒久的blob保存までは保証しない。重要な結果ファイル・コードは元Git blobを選択的に統合する。第三者のEEA/GBIF/EPA原アーカイブは本リポに再配布しない。",
      ""])
    Path("docs").mkdir(parents=True,exist_ok=True)
    Path("docs/PR_RESEARCH_REGISTRY_20261116.md").write_text("\n".join(doc),encoding="utf-8")
    print("GESUI_PR_REGISTRY_VERIFIED="+json.dumps(summary,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--token",default=os.getenv("GITHUB_TOKEN") or os.getenv("GH_TOKEN"))
    p.add_argument("--branch",default="analysis/gesui-all-evidence-one-plane-20261116")
    p.add_argument("--output",type=Path,default=Path("provenance"))
    main(p.parse_args())
