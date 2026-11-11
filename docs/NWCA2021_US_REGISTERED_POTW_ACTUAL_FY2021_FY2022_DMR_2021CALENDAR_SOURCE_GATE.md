# 全国2021年湿地調査 × 2021年実際の下水DMR報告：原本SOURCE-ONLY照合

**2026-11-11。分析の主線は全国・全球の生物群集であり、局所のシジミチョウ解析に戻さない。** この解析は、NWCA 2021の登録下水処理放流口近接と植物Shannon多様性の結果（等地点重み−0.12043、NWCA調査面積重み+0.18062）を**既に知った後の探索的な実排水証拠の確認**である。2021年の真の下水流入に直接対応していたとする前提は置かない。下水の効果を肯定または否定するモデルの事前登録ではない。

## 固定済みのsourceと応答を切り離す
- 元のEPA NWCA2021 **568**標準化湿地のSOURCE-ONLY SITE_ID + VISIT_NOを保持（実測水質と植物の対応が存在する調査区）。**植物多様性（H_ALL, H_NAT, H_ALIEN, FQAI）、水質（TN/TP）の値をこの段階で開かない。**
- 2026-11-08更新のEPA ECHO NPDES放流口ZIP（既に固定済み原本SHA256 `f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160`）。NAICS **221320**またはSIC **4952**、登録地理座標 `LATLONG_TYPE='Permitted Feature'`、台帳の許可有効日から2021年存在可能性を判定した厳格ソースに限定。2021年時点の実際の許可更新歴・全排水を完全に再現できない点を明記。
- 既にsource-onlyで登録放流口10km以内に**128湿地**（行政許可番号は**160**個）、州×生態地域の遠隔20km以上と比較可能な**240湿地・27州**（近103・遠137）を判定している。**これらの地点や距離を新たなDMR報告内容を見て選び直さない。**
- 今回新規に照合するのはEPA ECHO ICIS-NPDES **実際のDischarge Monitoring Report（DMR）** の **FY2021** アーカイブと **FY2022** アーカイブ。2021暦年全体を網羅するため、FY2021から2021年1–9月、FY2022から2021年10–12月を抽出。FYと暦年を混同しない。元の公表zip URL `https://echo.epa.gov/files/echodownloads/npdes_dmrs_fy2021.zip` と `https://echo.epa.gov/files/echodownloads/npdes_dmrs_fy2022.zip`。毎週更新される2026年抽出であるので、**抽出時の全ZIP SHA256**を保存。
- DMR原本の監視期間末日 `MONITORING_PERIOD_END_DATE` で選択。各permit `EXTERNAL_PERMIT_NMBR` の少なくとも1レコードがあることと、`PERM_FEATURE_TYPE_CODE=EXO`、採水地点 `MONITORING_LOCATION_CODE in {1,2,EG,Y}`、`DMR_VALUE_NMBR`が有限の数値、`NODI_CODE`空欄（no dischargeなどの未測定指標を実測と扱わない）を分ける。
- **TN：PARAMETER_CODE=00600**、**TP：00665** を最も厳格な既定値とする。それ以外のNO3/TKN/PO4形態を全窒素・全リンと合算しない。質量値Q1/Q2と濃度C1/C2/C3は別集計し、**単位不一致の実測濃度・質量をそのまま合算しない**。この段階では実測値の大きさをモデル化せず、適切なソース報告の*存在と欠損*のみ調べる。
- 公式元表の許可番号と放流口番号を照合。先行のsource-only JSONは湿地別に近接許可番号SHA256を持つが、実際の2021年報告が**同じ近接放流口feature ID**なのかは自明ではない。最初は**同じ行政許可番号にDMR報告がある**ことと、feature×許可番号が一致する厳格条件を**別々に集計**する。前者だけで後者の代用をしない。

## データ品質の判断を新規のTN/TPや植物の符号を見る前に固定
- この段階は**SOURCE-ONLY**。主たる機構比較を許可するには、元の128近接湿地のうち同一行政許可番号で2021暦年実測DMRが確認できる湿地**60以上**、対応する独立許可番号**30以上**、**10州以上**、そのうち**同じ登録feature（EXO）**の記録が確認できる湿地**40以上**、独立featureあるいは施設**20以上**を要求。実測TNおよびTP両方を持つsourceが**15施設以上**なければ、「観測された2021年N/P排水負荷で説明する」という解析は停止。少数の源を複数湿地や調査日で疑似反復させない。
- 何も記録されていない許可を**無放流**と判定しない。DMR報告は州による義務・規制の差があり、真の下水量の指標として一様ではない。水域に流れ着く**NHDPlusの方向付き接続**が未確認の地点を「下水曝露確認済み」には数えない。
- SOURCEゲートが通っても、2015年WWpの植物GBIF正相関や2021年NWCAの正の調査設計重みとの類似だけで機構を断定しない。NWCA植物H_ALLとの正式な栄養流入推定には**実測放流N/P（報告数値・単位・期間・流量）、湿地受水域、上下流水文接続**の独立検証が必要。
- 集計は全国の160原登録許可番号を母集団とし、選択は2021年の実排水結果に依存しない。公開アーカイブ原本はコミットせず、行政許可番号はSHA256化。過去の肯定結果も保持し、失敗または観察不足も残す。Draft PRのみ、main/chochoの原稿変更なし。

## 原データ
EPA ECHO ICIS NPDES DMR docs: https://echo.epa.gov/tools/data-downloads/icis-npdes-dmr-summary
EPA DMR national fiscal-year archives: https://echo.epa.gov/tools/data-downloads
Original outfall source frozen in: `results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json`.
