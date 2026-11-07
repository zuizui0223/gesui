# 非公式 HydroRIVERS ミラーへの暫定フォールバック

公式 HydroRIVERS Europe のZIPは GitHub Actions 実行環境で HTTP 403。取得制限を「生態学的な効果がない」と解釈せず、公開されている全球 HydroRIVERS v10 の第三者提供 ArcGIS FeatureServer（Forest Stewardship Council）を**暫定照合**に使用する。URL: https://maps.fsc.org/server/rest/services/Hosted/Optimized_Hyrdo/FeatureServer/0。

HydroWASTE は原本を SHA で固定し、植物記録を使わずに、WWTP・放流先・河川IDを選択。元河段 HYRIV_ID、上流を示す NEXT_DOWN の逆関係、下流の NEXT_DOWN を照会する。HTTPで取得した全レスポンスを SHA256 で監査し、利用サンプルを固定する。

**重要な限界**：ミラーは公式データの同一バイト原本と照合できていない。生物記録の取得前に、河段のつながり、サンプルの非重複と品質を確認する。処理水の体積は植物に届く N/P の実測質量ではない。GIS段階のみで植物群集の効果を結論しない。