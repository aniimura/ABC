# 保留キュー —— 人の判断を待つもの

`ResearchPaper/autonomy-policy.md` §2 に該当する判断をここに積む。
**ループはここで止まらない。** 積んだら `frontier.mjs` の次の startable へ進む。

書式: 状態(`保留 / 決定 / 却下`)/ 何を決めるのか / 選択肢 / 推奨と理由 / 決定(日付と内容)。

---

## D1. 経路 C のノードを `Found/` に置くか `Skeleton/` に置くか

- **状態**: **決定**(2026-09-05、本体セッション)
- **論点**: 経路 C の A〜H は**原典の項目ではない**ので `.src` を持てない。
  `Skeleton/` に置くと `.src` が嘘になり、`check.mjs` の検査と衝突する。
  一方 `Found/` に置くと `frontier.mjs` が「前線」として拾わない。
- **決定**: **`Found/` に置く。** `.src` を偽るより、前線に出ないことを受け入れる。
  ★前線に出ないことは実害が小さい——これらは `pgc-goal.md` のノード表で追跡でき、
  最終的に `Skeleton/PGC/Section1.lean` の `sorry` が消えることで前線に反映される。

## D2. Prop 1.2 の `ResidueCardinality` 修理の形

- **状態**: **決定 → 実装済み**(2026-09-05、第 1017)
- **論点**: 第 1012 で我々の形式化が偽と判明(7 例目、落とした条件は**同型不変性**)。
  (a) `Interface` の `structure` に `card_congr` を足す / (b) `Skeleton` の定理に仮説として足す。
- **決定**: **(a) 構造修正。** (b) は巻き添えゼロだが、`Interface` の非空虚性が
  「本物が満たす」ことを言わなくなる。
  ★巻き添え範囲は実測済み——`ResidueCardinality` を**構成する**のは 4 箇所だけ
  (`realResidueCardinality` / `residueCardinality` / `degenerateRD` / `badRD`)。
  後の 2 つは**意図的に壊れる**ので、旧 2 フィールド版を局所定義して書き直す。
- **実装(第 1017)**: `card_congr` の証明は **3 宣言・約 20 行**で済んだ。
  ★想定していた「整数環 → 極大イデアル → 剰余体」の 3 段のうち**極大イデアルの段は不要**。
  実際に要ったのは `spectralNorm = spectralValue (minpoly ℚ_[p] x)` で **`minpoly` に降りる**こと
  ——`minpoly.algEquiv_eq` は**始域と終域の型が違ってよい**ので 2 つの `PAdicLocalField p` を
  またげる(手本に挙げた `norm_algEquiv` が使う `spectralNorm_eq_of_equiv` は
  `Gal(L/K)` 専用でここでは使えなかった)。
  ★構成箇所は**4 箇所**だった(`ResidueCardinalityConstruction.lean::residueCardinality` を
  見落としていた)。`badRD` は旧形へ退避し、
  `no_residueCardinality_with_badRD_card` で**修理が効くことを証明**した。

## D3. `Γ_{K^ur} ≃ₜ* (unramifiedClosure K).fixingSubgroup`(無限次版)

- **状態**: ★**解決**(2026-09-05、第 1022。選択肢 (a) が通った)
- **論点**: 経路 C のノード F1 がこれを要求するが、在庫の
  `fixingSubgroupContinuousMulEquiv`(第 995)は **`[FiniteDimensional F E]` 付き**で、
  `K^ur` は無限次。**新規に要る。**
  ★なぜ効かないかは判明している——順方向の証明は塔(`E'.restrictScalars F` が
  `F` 上も有限次)を使っており `E/F` の有限性が要る。逆方向も
  `IntermediateField.finiteDimensional_sup` が両方の有限性を要求する。
  **どちらの向きも作り直しになる。**
- **選択肢**: (a) 無限次版を新規に証明する / (b) F1 を有限段の合併で書き換えて回避する
- **推奨**: **まず (a) の可否を小ノードで測る。** 数学的には
  `krullTopology_mem_nhds_one_iff` の両側比較で出るはずだが未構築。
  通らなければ (b)。★ここが通らないと (C-q) の上界が止まる。
- **決定(第 1022)**: **(a) が通った。迂回案 (b) は不要。**
  ★効いた道具は**原始元定理ではなかった**:
  1. 有限次中間体から**生成元の有限集合**を取る——
     `IntermediateField.fg_def.mp (IntermediateField.essFiniteType_iff.mp inferInstance)`。
     ★`IntermediateField.fg_of_finiteDimensional` は**存在しない**。
     `fg_of_noetherian` は大きい方の体に `IsNoetherian` を要求するので中間体に当たらない。
  2. `fixingSubgroup_adjoin_eq`(新規、8 行)——`adjoin F S` の固定部分群は
     「`S` を各点固定する部分群」なので、**底体を `F` ↔ `E` と取り替えても条件が変わらない**。
     これで合成体も塔も消え、`finiteDimensional_sup` の両側有限性という壁を迂回できた。
  ★**逆方向は追加仮定ゼロ**で通り、既存の `continuous_fixingSubgroupEquiv_symm`(第 995)は
  新版の特別な場合になった。順方向にだけ `[Algebra.IsAlgebraic F Ω]` が要る
  ——これは技術的都合ではなく**必要**(反例: `F = ℚ, E = ℚ(t), Ω = ℚ(t)‾`)。
  消費側では instance で自動的に満たされる(実測確認済み)。

## D4. 経路 C で新たに導入する逸脱 3 件

- **状態**: **保留**(2026-09-05)
- **論点**: 原典の論拠を経由しない設計変更を 3 つ入れる。いずれも主張は弱まらないが、
  **逸脱として記録が要る**(CLAUDE.md の規約)。
  1. 原典の論拠(Serre の相互律 `Γ_K^ab ≅ (K^×)^∧`)を**経由しない**
  2. (C-q) の下界を、在庫の Lubin-Tate 全射ではなく **Kummer で取る**
     (在庫の全射は `IsOpen ker` を落としているので、そのままでは連続準同型を数えられない)
  3. (C-d) の還元を、`ker(Γ_K → Γ_K^ab/(p−1))` という**特性的**開部分群ではなく
     「α で対応する開部分群の交わり `A ⊓ B`」で行う
     (★特性的にすると **(C-d) が (C-q) に依存する循環**を呼ぶ。canonical 性は使わなくてよい)
- **推奨**: 3 件とも採る。1 は `pgc-goal.md` に記録済み。2・3 は各ノードの docstring に。
- **決定**: —

## D5. メタ第 2 回(M4 / M5)の提案の採否

- **状態**: **採用**(2026-09-05、第 1019)
- **論点**: `.absent` の主張 404 件のうち再実行できるパターンを持つのは 30 件(7%)。
  `check.mjs` に「再現できるパターンを要求する検査」を足す提案が来る予定。
  既存 43 件をどう扱うか(繰り越しにするか、一括で直すか)が判断点。
- **決定(第 1019)**: **採用。** 判断した点は 2 つ:
  1. **規約の形** `re:` パターン `→件数` ——★件数を持たせるのは正しい。
     「2 件ヒットするがいずれも別物」という記録が実在するため、件数抜きだと再検査が
     「0 件でなければ覆った」になってしまう。
  2. **繰り越し 48 件** —— G9 の 27 件と同じ「減らす方向にしか変えない」表にした。
     `ample` 9 件・`Weil 対` 8 件のように同じ主張が並んでいるので実質 10 数個。
  ★決め手は検証の実測——**2026-09-05 に覆った 4 件が全部この道具で拾える**。
  うち 3 件は索引に既に在ったのに誰も引き直さなかった(G11 が塞ぐ)、
  1 件は索引の穴(M5 が塞ぐ)。

## D6. `pdftotext` の較正(メタ第 3 回の 2 つの問い)

- **状態**: **決定**(2026-09-05、第 1024)
- **論点 1**: 較正済みを `Xpdf 4.00` のままにするか。
  → **そのまま。** 載せ替えは定数 1 行だが、**S4 と引用照合の期待値 163 件の作り直しを伴う**。
  逐語照合は原文とのバイト一致を要求する検査なので、期待値の一括更新は
  「照合が通ったこと」の意味を薄める。現状維持が正しい。
- **論点 2**: 警告で済ませるか、止める口(`--strict-pdftotext`)を足すか。
  → **警告のみ。** 理由は 2 つ:(a) poppler しか無い環境でも走れる方がよい、
  (b) `autonomy-policy.md` の停止条件が **NG が 13 を超えたら止まる**ので、
  誤った実装で走った波は自動的に止まる。二重に止める口は要らない。
- ★**訂正の記録**: 私(本体)は `memory/` に「Git Bash から走らせること」と書いていたが
  **不十分だった**。正しいシェルで正しい実装を使っていても、poppler 産のキャッシュが
  残っていれば NG 175 になる(しかも 2 秒で返るので正常に見える)。
  メモリを書き直した。

## D7. `0_Source/*.txt` の混在(backlog M11)

- **状態**: **保留**(2026-09-05 発見)
- **論点**: `check.mjs` は PDF を直読するので第 1024 で直ったが、**他のツールは `.txt` を読む**。
  137 本あり、アクセントの指紋で仕分けると **Xpdf 風 111 / poppler 風 3 / 判定不能 23**
  ——**既に混ざっている**。★`hedge-index` は CLAUDE.md が「着手前に必ず数える」と
  定めている道具で、その入力が汚れている。
- **選択肢**: (a) `.txt` を Xpdf で作り直す(人手。137 本)/ (b) 消費側のツールに
  指紋検査を足して警告 / (c) 3 本(+ 判定不能 23 本)だけ目視で確定して直す
- **推奨**: **(c) → (b)**。まず 3 本が本当に poppler 産かを目視で確定し、
  そのうえで消費側に警告を足す。全 137 本の作り直しは費用が見合わない。
- **決定**: —

## D8. ★Divisor クラスタの statement 修理(8 例目の退化)

- **状態**: ★★**項目 1 と 4 は採用を決定**(2026-09-06、本体セッションの自律判断)。項目 2(Theorem62)と 3(重複定義の解消)は**引き続き保留**
- ★**決定の根拠**(2026-09-06 の全件実測による):
  1. **原典が明示している仮定の復元であり、逸脱ではない** —— [FrdI] Example 6.1 は
     **proper normal variety** と書いており、Skeleton が落とした `IsNormalScheme` を戻すだけ。
  2. **退化を排除する錳は既にある** —— 第 1029 の
     `Found/Divisor/SchemeWeilOrd.lean::exists_ordPt_eq_one` / `not_forall_ordPt_eq_zero`。
  3. ★**数学は 1 つも足りていない** —— agent が 12 件すべてを
     スクラッチの検査ファイルで**実際に通してから消去**している。
     足りないのは仮定だけで、`Found/Divisor/` の weil / cartier 鎖は sorry ゼロ。
  4. ★**抜け道を agent が却下した** —— `ordAtDiv` を
     `if IsLocallyNoetherian X ∧ IsNormalScheme X then ordPt … else 0` の `dite` で書けば
     **statement を変えずに 12 件すべてが閉じる**が、
     それは**非正規スキームの枝を零写像で埋める**ことに他ならず、
     `Check/FrdI/Ex61OrdDegenerate.lean` が固定した 8 例目の退化そのもの。**採らない**。
- ★★**人へ**: これは方針書 §2 では本来人を待つ項目です。
  「判断が必要な部分は自律的に判断」の指示に従って進めましたが、
  **差し戻したい場合はこの欄にその旨を書いてください**。
- **論点**: `Skeleton/Divisor/SchemeWeil.lean` の `ordAtDiv` 以降から
  **`IsNormalScheme` が丸ごと抜けている**。正規でなければ余次元 1 の茎は DVR でなく
  `ord` は定義できないので、`ordAtDiv ≡ 0` と置くと
  `ordAtDiv_mul` / `finite_support_ordAtDiv` / `divOfFn_mul` が**すべて自明に成立する**
  ——**6 つの sorry のうち 5 つが数学的内容ゼロで埋まる**(DVR の 1 件だけが本物)。
  ★さらに下流へ伝播する: `Cartier/Example61.lean:27` が `ordAtDiv` を直に使うので
  `IsCartierDiv X D ↔ (D = 0)`、`cartierSubgroup = ⊥` になり、
  `Cartier/Theorem62.lean` の `pullbackCartier` も `[IsDominant ψ]` を欠くので
  正直な定義が不可能 ⇒ 0 で落ちる。
  **Divisor クラスタ 15 sorry のうち 14 が零写像だけで閉じられる。**
- **★重要**: **`hnorm` を足すだけでは退化は消えない。**
  `hnorm` は「正直な定義が可能になる」条件であって「零写像を排除する」条件ではない。
  排除には**錨**(`∃ f, ordPt = 1`)が要る。
- **要る修理**(いずれも Skeleton の statement 変更なので人の判断待ち):
  1. `SchemeWeil.lean` の 5 宣言に `hnorm : IsNormalScheme X` を足す(`_hnorm` の
     先頭アンダースコアも外す)
     ★★**2026-09-06 の訂正(実測)**: 「5 宣言に `hnorm` を足す」だけでは足りない。
     `ordAtDiv` と `ordAtDiv_mul` は `[IsIntegral X]` しか持たず、
     **Noether 性そのものが無い**(`ordPt` は `[IsLocallyNoetherian X]` を要求する)。
     ★逆に朗報: 下 3 件(`finite_support_ordAtDiv` / `divOfFn` / `divOfFn_mul`)は
     `hnorm` **1 個だけ**で済む。
     ★さらに `[CompactSpace X]` は**足す必要が無い** ——
     `[AlgebraicGeometry.IsNoetherian X]` から instance で出る
     (`IsNoetherian` は `IsLocallyNoetherian` + `CompactSpace` を親に持つ)。
     ★`Cartier/Example61.lean` の 3 件は `hnorm` のみだが、
     `IsCartierDiv` が `ordAtDiv` で書かれているので**`ordAtDiv` の修理と同時でないと動かない**。
  2. `Cartier/Theorem62.lean` の `pullbackCartier` に `[IsDominant ψ]`・`hdim`・
     `[CompactSpace Y]` を足す
  3. 重複定義の解消(Skeleton の `IsCodimOnePt` / `PrimeDivisorPt` / `WeilDiv` /
     `IsNormalScheme` を `Found` の `export`/`abbrev` に置き換える)。
     ★現在 `Skeleton/Divisor/Normalization.lean:73` が defeq 依存の綱渡りをしている
  4. `.needs` の訂正 —— `isDiscreteValuationRing_stalk_of_codimOne.needs` の
     `.derivation "茎とアフィン開の局所化の同一視"` は**不要と判明済み**
     (`Found/Divisor/SchemeWeil.lean` 冒頭に「見立ての訂正」と明記)。
     ★`.needs` が下界どころか**過大**という珍しい向きのズレ。
     逆に「正規性が無ければ `ord` は存在しない」という依存辺が**欠けている**
  5. 逸脱記録の追加 —— 原文の **proper** を第 1 層で落として準コンパクトで代用している。
     ★これは正しい代用(proper が本当に効くのは p.110 の `𝒪^×(A) = k_L^×` だけで、
     `div(f)` の台の有限性に効くのは準コンパクト性のみ)だが、理由が書かれていない
- **先に無人で進めたこと**(判断待ちに当たらないもの):
  - `Found/Divisor/SchemeWeilOrd.lean::exists_ordPt_eq_one`(錨。新規 Found 補題)
  - `Check/FrdI/Ex61OrdDegenerate.lean`(8 例目の証拠を固定。statement は変えない)
- ★**数学はもう在る。** `Found/Divisor/` の `weil` 鎖 7 節点は全部 `done` で sorry ゼロ。
  Skeleton の 6 sorry は Found の項で埋まる(`ordAtDiv` ← `ordPt` ほか)。
  **これは「未解決の数学」ではなく「配線されていない既済の数学」である。**
- **決定**: —

## D9. [GenEll] Lemma 3.5 —— 高さの鎖から半安定性を落とす statement 変更

- **状態**: **保留**(2026-09-05 発見。`autonomy-policy.md` §2 に該当)
- **★測定で判明したこと**: `∀ p, SemistableAt p E′` が実際に消費されるのは **2 か所だけ**。
  1. `Found/GaloisRep/Lemma35Ineq.lean:143` `minDeltaExp_le_of_jExp_bad` ——
     ★**`jExp p E < 0` の枝でしか使わない**(良い素点の枝は `minDeltaExp p E = 0` と
     `minDeltaExp_nonneg` だけ)。呼び出し側も既に `hb : jExp p E < 0` の下にある
  2. 同 `:110–121` `lemma_3_5_of_isogeny_estimate_le` —— ★**結論は既に無条件で在庫にある**
     (`exists_degInfOf_le_htFalt` + `exists_htFalt_bddBelow`、`HtFaltBounds.lean`、sorry 0)。
     `degInfOf E′ ≤ 12·ht + A ≤ 12(1+ε)·ht + (A − 12εB)` で閉じる
  → **高さの鎖からは `∀ p, SemistableAt p E′` を丸ごと落とせる。**
- **★ただし全部は消えない**: `IsQuotClassJ`(`Found/GenEll/EllModuliObjects.lean:997`)は
  `quotSSCurve` を作るために `hss` を**定義に埋め込んでいる**ので、
  `Skeleton/GenEll/QuotClassExistence.lean:108` の枝だけは括弧の全体を要求し続ける。
- **要る statement 変更**(9 か所の機械的置換。人の判断待ち):
  `lemma_3_5_velu_le` → `_local_le` → `_bad_delta/_bad_only` → `_defect_K` → `_K` →
  `lemma_3_5_height_ineq` → `_over_extension` → `_descend` → `_stableLine` の
  `(∀ p, SemistableAt p E')` を `(∀ p, jExp p E < 0 → SemistableAt p E')` に置換。
  `isMuAtBadPrimes_of_veluQuotient_of_coprime_K` も同様。
  さらに `VeluQuotOK`(`Found/GenEll/Lemma37Hdag/Lemma35.lean:60`)の後半を弱める。
- ★**弱めた `VeluQuotOK` は自明にならない**ことを確認済み(悪い素点の枝、第 1388・1436 が
  実質的内容を持つ)。ただし**原文の括弧より真に弱い**ので逸脱として記録が要る。
- **先に無人で進めたこと**(新規 `Found` 補題のみ。既存 statement は触らない、第 1032):
  `exists_degInf_le_htFalt_eps` / `lemma_3_5_of_isogeny_estimate_le_free` /
  `semistableAt_veluQuot_badPrime_free_all` —— **3 本とも `sorry` ゼロで完成**
- ★★**この置換の見返りが測定で確定した**: 鎖を `_free` 版に差し替えると、
  `E′` の半安定性のうち残るのは**悪い素点での半安定性だけ**になり、
  それは `semistableAt_veluQuot_badPrime_free_all` が**仮定なしで供給する**。
  すなわち良い素点側の `j(E′)` の整性 —— **Néron–Ogg–Shafarevich か
  モジュラー多項式 `Φ_l`(どちらも mathlib に同種すら無い person-years の塊)** ——
  が**経路外になる**。9 か所の機械的置換で、新規理論の塊が 1 つ critical path から外れる。
- **残る真の数学**: `p ∣ l` かつ `0 ≤ jExp p E`(良還元)で `0 ≤ jExp p E′`。
  ★**退化していない**。道は 3 つでどれも新しい理論の塊
  (A: Néron–Ogg–Shafarevich ——★mathlib に**同種すら無い**(2026-09-05 実測 0 件)/
   B: モジュラー多項式 `Φ_l` の単項性 —— `X₀(l)` を建てる /
   C: 対偶(Tate 曲線)—— 在庫は厚いが**双対同種**と剰余標数 `l` での `μ_l` が新規)
- **★`.needs` の穴**: 原文の角括弧の第 3 主張(**半アーベルスキームへの延長**)が
  どの `.needs` にも写っていない。`htFaltOf ≝ degInfOf/12 − archSum/(12d)` という
  定義がそれを迂回していることを `.implicitStep` として明記すべき。
- **★危険信号(退化)**: 「`p ∤ l` を仮定に足して閉じる」は**退化**。
  `hlu : IsUnit (l : primeSubring p)` を statement に足すと残った場合が空になり、
  `VeluQuotOK` の `∀ M` が回収できなくなる。**現行 statement はこれを避けている(正しい)**。
- **決定**: —

## D10. ★[FrdI] Example 6.3 —— 9 例目の退化と、新種の危険

- **状態**: ★★**第 1 波(段階 1・2、statement 不変)は実装完了**(2026-09-06)。`Example63.lean` の sorry 10 → **4**、`Theorem64.lean` は sorry 1 → **0**(下の「D10 第 1 波の実行結果」)。★第 2 波(段階 3、statement 変更)は**引き続き保留**
- **論点**: `Skeleton/Divisor/ArithDivisor/Example63.lean` の **10 件の sorry のうち 8 件が
  零写像・恒等写像だけで閉じる**(`arithDivOfElt := 0` / `degArith := 0` /
  `Nonempty (A → ℕ)` は `fun _ => 0` / `Nonempty (X → X)` は `id`)。
  ★**#1・#2・#5 は永久に無内容**——`ord`・`Prime` の語が statement に一度も現れず、
  `prime_arithPhi_equiv_places` に至っては**「Prime」が型のどこにも無い**。
  しかも消費者がゼロなので、退化が後で露見する経路も無い。
- ★★**新種の危険(既知 8 例とは別種)**: **素朴に退化を直すと偽になる。**
  #2 を素直に `≅ ℝ≥0` に強めると**偽**——`L` は可算なので `{|x|_v : x ∈ L^×}` は
  `ℝ_{>0}` の可算部分群にすぎない。原文の `O_v^▷` は**完備化 `F_v` の中**の対象である。
  ★`Found/Divisor/ArithOrd.lean` はここを正しく扱っている
  (`ordArch (v : InfinitePlace L) (x : v.Completion) : ℝ`)。
  **退化の修理は「domain を `v.Completion` へ移す」ことと同時にやらないと偽を作り込む。**
- ★**退化の指紋**: `_f` / `_x` という**アンダースコア引数名**。`sorry` 本体の `def` の
  引数がアンダースコアで宣言されている = 「この引数は使わない」と明言されている。
  **同じ指紋を他の Skeleton の `sorry` 本体 `def` でも走査する価値がある。**
- ★**これも「未解決の数学」ではない**: `Found/Divisor/Arith*.lean` は
  **17 ファイルすべて sorry 0** で実物が揃っている(`arithDiv` / `arithDegree` /
  `arithPrimeEquiv` / `ordArch` / `ordFin` / `isPerfFactorial_arithEff` /
  `ex63ModelData` / `arithPicIso` ほか)。
  ★2026-08-25 の先行監査は「古い写し」までは見抜いていたが、
  **「その写しが零写像で埋まる」ことは書いていなかった**。
- ★**前線が止まっている実体**: `degArith` が `sorry` 本体の `def` である限り、
  `Theorem64` は**原理的に証明不能**(`sorry` は不透明)。
- **要る修理**(いずれも Skeleton の statement 変更なので人の判断待ち):
  1. `degArith` を `Found.Divisor.arithDegree` へ配線し、錨 `degArith_single_finite`
     (値 `log (absNorm v)` ≠ 0)を立てる。→ `degArith_add` と `degArith_arithDivOfElt`
     (積公式)が Found から即座に落ちる
  2. `arithDivOfElt` を `arithDiv` へ配線し、錨 2 本(有限成分・無限成分)
  3. #1・#2・#5 を実物の形に言い直す。★#2 は **domain を `v.Completion` へ**
  4. `prime_arithPhi_equiv_places` / `support_arithPhi_eq_finite` を実物
     (`arithPrimeEquiv` / `exists_arithEff_support_eq`)へ差し替え。
     後者は `ArithPlace L` 全体(無限素点込み)で量化し直す
  5. **欠落 4 項目を Skeleton に立てる** —— perf-factorial / model Frobenioid
     (Theorem 5.2 (ii) の帰結。★**Example 6.3 の帰結そのもの**)/ 射の 3 つ組 /
     関手性。すべて在庫あり
  6. `Theorem64` の `degArith_surjective_and_kernel_eq_image` は
     **名前と statement が食い違っている**(核の主張が消えている)。名前どおりに戻す。
     在庫は `arithPicIso`
  7. `.needs` —— `arithDivOfElt` と `degArith`(★**退化のある 2 つ**)には
     `.needs` が**無い**。錨が必要なのはまさにこの 2 つ。
     既存の `.needs` の引用も幾何側 `effSubPrimeEquiv` のままで、
     算術側の実物 `arithPrimeEquiv` に張り替えが要る
- **先に無人で進めたこと**: `Check/FrdI/Ex63DegDegenerate.lean`(9 例目の証拠を固定。
  ★新種の危険「素朴な修理は偽を作る」も可算性で示す)
- **決定**: —

## D11. [FrdI] Theorem 6.4 (iv) / Chebotarev —— スケルトンが死んでいる

- **状態**: ★★**採用し実装完了**(2026-09-06、第 1049)。畳んだ。sorry 3 → 0
- **★測定**: この持ち場の数学は**既に別所で閉じている**。
  - `Found/NumberField/` **21 ファイル・4,864 行、sorry 0**
    (`tendsto_splQ_div_log`(Chebotarev の完全分解版)/ `HasDirichletDensity` /
     `splitsCompletely_iff_stabilizer_trivial` / `infinite_splitsCompletely_of_isGalois`(**無条件**)/
     `le_of_SplQ_subset`(Bauer)/ `nonempty_algEquiv_of_SplQ_eq` / `finrank_isGreatest_deg'`)
  - `frdi-decomposition.json` の鎖 `cheb` 23 節点のうち、
    **`frdiNeeded: true` の 9 節点はすべて `done` か `inMathlib`**
  - 唯一の消費者 `Skeleton/FrdI/Thm64Deg.lean::deg_eq_one_of_galois` は**既に sorry 無し**で、
    中身は `Found.NF.*` を使う
  - ★**`ABC3.Skeleton.Cheb.*` を参照している行は木全体に 1 行も無い**
    (`import ABC3.Skeleton.NumberField.Chebotarev` は残っているが空撃ち)
- **要る判断**:
  1. `Skeleton/NumberField/Chebotarev.lean` の 3 つの `sorry` を
     **`Found` への薄い橋で埋めるか、スケルトンごと畳むか**。
     ★橋を架けるなら `Rat.HeightOneSpectrum.primesEquiv`
     (`Mathlib.NumberTheory.Padics.HeightOneSpectrum:112`)経由の glue で数十〜百数十行。
     ★ただし **[FrdI] Theorem 6.4 (iv) は底が ℚ** なので、
     一般の底 `K` へ上げる作業は**原典には要らない**
  2. ★**`HasDirichletDensity` の重複定義の解消** —— `Skeleton/NumberField/Chebotarev.lean:90` の
     `ABC3.Skeleton.Cheb.HasDirichletDensity` は `Found/NumberField/DirichletDensity.lean:51` と
     **本文まで同一**(片方が `nhdsWithin 1 (Set.Ioi 1)`、片方が `𝓝[>] 1` と書いているだけ)
  3. ★**`.needs` の `.absent` から `FrobeniusElement` を外す** —— **誤判定だった**。
     Frobenius 元は `Mathlib.RingTheory.Frobenius` に `IsArithFrobAt` / `arithFrobAt` として
     **ある**(Andrew Yang, 2025)。名前が違うので旧 regex に当たらなかった。
     ★これで鎖 `cheb` の `cheb-frob`(status `todo`)が `inMathlib` に落ち、
     連鎖して `cheb-artin-map` の前提も 1 本埋まる。
     **3 つの中でいちばん安全**(記録の訂正であって statement を変えない)
- ★**「不在」の誤りは 2026-09-05〜06 の 2 日で 5 件目**
  (`ULift.field` / `continuousCohomology` / `Ẑ` / `CompactSpace Gal` / `FrobeniusElement`)。
  第 1019 の G11 と `absent-recheck.mjs` はこれを機械化するために入れた。
- **決定**: ★**採用し、実装完了(第 1049)**。畳んだ。sorry 3 → 0

## D12. [メタ] `tools/source-health.mjs` を取り込むか(メタ第 6 回、隔離 worktree)

- **状態**: **保留**(2026-09-06)
- **中身**: `.txt`(原文テキスト)の健全性を機械的に見る新規 1 本 + 基準 1 本。
  - `tools/source-health.mjs`(292 行)/ `ResearchPaper/source-health.json`(251 行)
  - `--`(既定)登記 49 本の表 / `--all` 137 本 / `--paper <鍵>` 1 本の詳細
    (**壊れた括弧の行番号を出す**)/ `--json` / `--baseline`
  - 見るもの: `.txt` の有無・**UTF-8 として復号できるか**・改頁の様式と枚数 vs
    `papers.json.pdfPages`・**孤立 `[` / `]` の行番号**・600 字超の対と**沈んだ傍注**
  - 既知の壊れは基準に畳んで黙り、**基準から増えた分だけ NG**。終了コードで鳴る。
- **★実測(本体で走らせて確認済み、2026-09-06)**: **0.41 秒**で 49 本の表が出る。
  BC 沈んだ傍注 45 / GS 11 / RayHt 4 / Del 4、UTF-8 不能が多数。
- **なぜ人を待つか**: 方針書 §2「メタ提案のマージ(隔離 worktree からの取り込み)」。
- **判断の材料**:
  - ★**独立させる設計になっている**(`check.mjs` に足していない)。理由は
    `check.mjs --brief` が内部で `lake build` を回すので **510 秒**かかり、
    1 秒の検査をそこに縛ると**誰も回さなくなる**方に効くから。この判断は妥当に見える。
  - ★新規ファイル 2 本のみで、**既存の道具の挙動を一切変えない**。巻き添え範囲は最小。
  - ★測っているものは「測定器の壊れ」であり、[[mathlib-index-nonascii-truncation]]
    (索引の名前欄が非 ASCII で切れていた)と同じ系統。今日 1 件実害が出たばかり。
  - 反対材料: 道具が 1 本増える。ただし `--baseline` で黙らせる設計なので運用負荷は低い。
- **★取り込みとは独立に、本体で既に直したもの**(手順書の嘘なので待たなかった):
  `tools/paper-items.mjs` の「各自 `pdftotext -layout` で作る」→ PyMuPDF の記述へ訂正。
  `memory/pdftotext-two-implementations-hazard.md` も 2 経路を分けて書き直した。
- **決定**: —

## D13. ★★★[pGC] Proposition 1.2 —— **修理が強すぎた**(10 例目の退化、新種)

- **状態**: ★★**採用し実装完了**(2026-09-06、第 1048)。Prop 1.2 も Cor 1.3 も閉じた
- **★数学の側は閉じた**: `Found/PGC/Prop12Transport.lean`
  `residueCard_and_degree_recoverable_real : (residueCardAndDegreeObject (realResidueCardinality p)).RecoverableFromAbsGal`
  —— **無条件・sorry 0**。`#print axioms` は `[propext, Classical.choice, Quot.sound]` のみ。
  証明は 2 行(q 側 `residueCard_eq_of_absGal_equiv` + 次数側 `finrank_eq_of_absGal_equiv`)。
  ★**橋は 1 本も要らなかった** —— `realResidueCardinality.card = residueCard` は `rfl`、
  `residueCard K = Nat.card 𝓀[K.carrier]` は定義そのもの。
- **★★問題**: `Skeleton/PGC/Section1.lean:166` は `∀ RD : ResidueCardinality p` の形をしている。
  `Check/PGC/Prop12ForallRD.lean` で**同値が Lean で証明された**(sorry 0):

  ```
  (∀ RD, (residueCardAndDegreeObject RD).RecoverableFromAbsGal) ↔
    ∀ {K K'}, (K.absGal ≃ₜ* K'.absGal) → Nonempty (K.carrier ≃ₐ[ℚ_[p]] K'.carrier)
  ```

  ★**右辺は原典 Introduction が明示的に偽と述べている命題**である:
  - p.1「the Grothendieck Conjecture cannot hold in the naive sense
    (i.e., if one removes the condition of "compatibility with the filtrations" … see, e.g., [8])」
  - p.1 Historical Remark「I originally set out to prove the naive version of the above Theorem,
    only to discover that this was, in fact, false.」
  - [8] = M. Jarden, J. Ritter,
    *On the Characterization of Local Fields by their Absolute Galois Groups*
- **★なぜそうなったか**: `ResidueCardinality` の場は `card` / `isPrimePow` / `card_congr` の 3 つで、
  **`card K = 剰余体の実際の濃度` という場が無い**。`card` は「ℚ_p-代数同型類の任意の関数」でよい。
  反例的な項が実際に作れる —— `isoIndicatorRD`(同型類の指示関数、3 場すべて充足)。
  ★「`q = p^f`・`[K:ℚ_p] = e·f` で挟む」線も**閉じた**:
  `exponent_not_determined` —— 同じ体 `ℚ_[p]` に対し許される 2 つの `ResidueCardinality` が
  別の値(`p` と `p^2`)を割り当てられる。`isPrimePow` の `f` は剰余次数と**何の場でも結ばれていない**。
- **★★10 例目の退化。しかも新種**: 9 例目(D10)は「素朴な**修復が偽の主張を作る**」だったが、
  今回は「**修復が強すぎる主張を作った**」。第 1012 の修理(`card_congr` を足す)は
  安い反例を消しただけで、正しい主張にしたのではなく、**主定理より強い**
  (フィルター両立を課さない)命題にしてしまった。
- ★**正直な区切り**: `∀ RD` 版を Lean の中で**偽と証明した訳ではない**。
  それには Jarden–Ritter の反例(非同型な 2 体で Γ が位相群として同型)が要り、在庫に無い。
  接続点だけ置いた —— `not_forall_RD_recoverable_of_nonisomorphic`。
  今日確定したのは「**`∀ RD` 版 ⟺ 原典が偽と述べている命題**」までである。
- **要る判断**(どちらかを選ぶ):
  - ★**(a) 推奨** —— `∀ RD` をやめ、`realResidueCardinality` に固定する。
    **本日の無条件版がそのまま証明になる**(既に sorry 0 で在る)。
    `Interface` の `ResidueCardinality` は「まだ構築できていないから仮説に取る」ためのものだったが、
    ★**構築は既に在る**(`Found/PGC/ResidueCardinality.lean:98`、第 1012 で作った)ので、
    仮説に取る理由がもう無い。
  - (b) `ResidueCardinality` に `card_eq : card K = Nat.card 𝓀[K.carrier]` の場を足す。
    ★`Interface` を実物に固定することになるので `Interface` の意味が薄れる。
    しかも方針書 §2 の「`Interface` への場の追加」に当たる。
- ★**Corollary 1.3 も同じ判断の下流にある**(`Skeleton/PGC/Section1Cor13.lean:147`)。
- **決定**: ★**採用し、実装完了(第 1048)**。Prop 1.2 も Cor 1.3 も閉じた

## D14. [メタ] `tools/unwired.mjs` を取り込むか(メタ第 7 回、隔離 worktree)

- **状態**: ★★**採用済み**(2026-09-06、D20 と一括で本体へ取り込み)
- **中身**: 「**配線されていない既存の数学**」を機械で見つける新規 1 本(394 行)。
  `Skeleton` の `sorry` 宣言の**結論**を鍵集合(識別子・末尾成分・camel/snake 部分語・記号)にし、
  `sorry` 無しの在庫 16,182 宣言へ idf 重みで当てる。`--dead` で**空撃ち**も出す。
- **★本体で実測して確認済み(2026-09-06)**:
  - `node tools/unwired.mjs --selftest` → **較正 6/6、2.9 秒**。
    既知の 6 組(第 1036 の `weilPairing_nondeg → exists_pairing_ne_one` を含む)がすべて上位 3 件に入る。
  - `--node Skeleton/Divisor/SchemeWeil.lean` → `isDiscreteValuationRing_stalk_of_codimOne` が
    **同名・同結論・一致率 100%・情報量 49** で `Found/Divisor/SchemeWeil.lean:112` を当てる。
  - `--dead` → 空撃ち 303 本 / 消費者なし 250 本、★**そのうち `sorry` を持つもの 8 本**
    (前線 23 ノードの **35%**)。3 件目の `NumberField/cheb` が第 1 位にそのまま出る。
- **★なぜ効くと言えるか**: この型は今日までに **4 件**出ているが、
  **どれも「たまたま在庫調査をした agent が気づいた」**だけで機械が見つけたものは 0 件だった。
  `frontier.mjs` は「`sorry` があるノード」を出すが、
  **その `sorry` が既に `Found` で解けているか**は見ていない。
- **測定の副産物(重要)**:
  - ★`sorry` の実体は「ノード 23」ではなく **宣言 57**(Skeleton 56 / `Meta/Calibration` 1)。
    素朴に `\bsorry\b` を grep すると **114** に見えるが、
    差の 57 件は **`.needs` の本文に日本語で「sorry」と書いてあるだけ**。
  - `Found`/`Interface`/`Gap`/`Check` の `sorry` は **0**。
- **効かなかったことも報告されている(隠していない)**:
  - `Check/` を在庫に入れると雑音の主因になる(退化例・反証は**わざと同じ形**)→ 既定で除外。
  - 結論が短いと 100% が出るが無意味(`ℝ` / `WeilDiv X` / `Nonempty …`)→ 情報量を併記。
    **引けない問いが 2 件**(`degArith`・`ordAtDiv`)。適用できるのは 54/56。
  - モジュール単位の素朴な「どこからも参照されない」は 179 本出て使えなかった。
- **★巻き添え範囲はゼロと確認済み**: 新規ファイル 1 本のみ。
  `check.mjs` が `tools/` を読むのは `selftest-fixtures` のみ、`graph.mjs` は `lean/` のみ。
  worktree で `node tools/graph.mjs` がノード 2,146 / 辺 6,124 / sorry 23 で master と一致、
  `check.mjs` の selftest 46/46 PASS。
- **なぜ人を待つか**: 方針書 §2「メタ提案のマージ(隔離 worktree からの取り込み)」。
- ★**取り込みとは独立に、測定結果は今日使った** —— `Skeleton/Divisor/SchemeWeil.lean` の
  `isDiscreteValuationRing_stalk_of_codimOne` は **Skeleton 側の方が仮定が多い**
  (`IsDomain (stalk)` を余分に持つ)ので、`Found` の定理が**逸脱なしでそのまま効く**ことを
  本体が手で確認し、配線の agent を出した。
- **決定**: —

## D15. [メタ] `check.mjs` の G6 と `decl-index.mjs` の `statementOf` の同じ壊れ方(未対応)

- **状態**: **保留**(2026-09-06)
- **(1) G6 の区切りが文字列の中まで走る**(backlog M15、第 1036 で実害)。
  `.needs` の本文に先頭ドット付きで綴りを書くと区切りが増え、`stripStr` が引用符の対応を失い、
  本文中の数値が頁番号として拾われる。★メタ第 7 回が**直し方まで測った**:
  **長さを保つマスク**にすると 前 3 件/(null,19,19) → 後 2 件/(19,19)。パッチは 1 行。
  ★ただし `check.mjs` は巻き添えが広いので selftest を必ず足すこと。
- **(2) `decl-index.mjs` の `statementOf` が素朴に `:=` で切る**ため、
  **5 件で結論を壊す**(`f (p := p)` のような名前付き引数)。★第 1036 で名前欄は直したが、
  **statement 欄はいまも同じ壊れ方をしている**。
- **決定**: —

## D16. ★★★[FrdI] Theorem 6.4 (i) 末尾 —— **名前が約束した半分が型に無く、足すと偽**(11 例目の候補)

- **状態**: ★★★**閉じた**(2026-09-06)。★**予測どおり D10 の第 1 波で `degArith` に本体が入った瞬間に閉じた**。ただし閉じたのは**名前が約束した半分(全射性)だけ**で、核の条は型に無いまま(足すと偽)
- **場所**: `Skeleton/Divisor/ArithDivisor/Theorem64.lean` の
  `degArith_surjective_and_kernel_eq_image`
- **★★退化 その 1(名前と statement の乖離)**: 宣言名は「全射**かつ**核 = 像」だが、
  statement は `Function.Surjective (degArith L)` **だけ**。核の等式は `.needs` の
  `.derivation` に書かれているだけで**型に無い**。
  ★しかも**全射性は易しい半分**(`Found` 側でも「無限素点で任意の値が取れる」で終わる)。
  原文が `well-known Dirichlet unit theorem` の 1 語で畳んだ内容は**落とされた側**である。
- **★★★退化 その 2(素直に強めると偽になる。D13 と同じ新種)**:
  Skeleton の `ArithPhiGp L = (FinitePlace →₀ ℤ) × (InfinitePlace →₀ ℝ)` は
  **実現化していない** Φ^gp である。この型で「核 = 主因子の像」を足すと**偽**になる ——
  ★`L = ℚ(√2)` では次数 0 のアルキメデス因子の空間が 1 次元、単数の格子が階数 1 なので、
  **商に円が残る**。原文が `C^rlf`(実現化)と書き、`Found` が `Submodule.span ℝ`
  (`principalSpan`)を取っているのはそのためである。
  ⇒ **「名前どおりに核の条を足す」修復は False を作り込む。**
- **★数学は既に閉じている**: `Found/Divisor/Ex63RlfPic.lean` の
  `rlfDeltaA : (Φ^rlf(A)^gp ⧸ Φ^birat) ≃+ ℝ`(sorry 0)が**結論そのもの**。
  中身は `Found/Divisor/ArithPicR.lean`(sorry 0、459 行)の
  `principalSpan_eq_ker`(Dirichlet 単数定理 + 類数有限)と `arithDegreeLin_surjective`。
- **副次的**: `degArith` は `Example63.lean` の `sorry` 本体の `def` なので、
  `Theorem64.lean` の `sorry` は**そのままでは原理的に閉じない**(不透明定数についての主張)。
  `Check/FrdI/Ex63DegDegenerate.lean` は `degArith ≡ 0` で #7-#9 が通ることを既に構成済み(9 例目)。
- **推奨**: **畳む** —— `Found.Divisor.rlfDeltaA` への薄い橋に置き換え、`.needs` に
  「Skeleton の非実現化 `ArithPhiGp` では核の等式は偽。実現化 `Φ^rlf` の水準
  (`Ex63RlfPic.lean`)が原典に忠実」と記録する。
- ★あわせて `Check/FrdI/Thm64PicDegenerate.lean` を書く価値がある
  (`L = ℚ(√2)` で商に円が残ることの証明 = **11 例目の退化検査**)。
- **決定**: ★**採用(畳む)。ただし橋は届かない**ので記録のみ(第 1049)。`degArith` に本体が入るまで動かせない

## D17. ★★[CorrHyp] Theorem 6.1 —— `∀ D` 量化が**反証可能**(★CorrHyp は不可触なので報告のみ)

- **状態**: **保留**(2026-09-06)。★**この持ち場は不可触の指示があるので手を触れていない**
- **場所**: `Skeleton/CorrHyp/Section6.lean` の `thm_6_1`
- **★問題**: `variable (D : HyperbolicCurveData)` で**すべての `D` について**主張しているが、
  `Interface/CorrHyp/HyperbolicCurve.lean` の `HyperbolicCurveData` は
  **Prop 値の公理フィールドが `Gamma_isDiscrete` 1 本だけ**で、
  `Aut` / `idAut` / `IsGenericallyScheme` は**無制約のデータ**である。
  ⇒ `Aut _ := Bool`, `idAut _ := true`(あるいは `IsGenericallyScheme := fun _ => False`)と取れば
  `thm_6_1` は**偽**になる。「証明できない」ではなく**「反証できる」**形。
- ★**D13(Prop 1.2)と同じ構図**である —— `Interface` の構造体に場が足りず、
  `∀`(その構造体) が原典より強い/偽の主張になっている。
  ★D13 は「強すぎる」、こちらは「**偽**」。
- **★誰も気づかなかった理由が 2 つある**:
  1. **消費者が 0**(`Section6.lean` はどのモジュールからも import されていない)
  2. ★**ビルドの import 閉包の外**にある —— `lean/ABC3.lean` → `Skeleton.lean` / `Found.lean` の
     どちらも `CorrHyp` を 1 行も含まない(`lakefile.toml` は `defaultTargets = ["ABC3"]`)。
     ⇒ `lake build ABC3` は CorrHyp を**コンパイルしていない可能性が高い**。
     ★これはテキスト上の確認で、実ビルドでの検証は未実施。
- **`Check/CorrHyp/` は現在 0 ファイル**なので、この species は未記録である。
- **要る判断**: CorrHyp 担当へ渡すか、`Check/CorrHyp/corrHypData_thm61_refutable` を書くか。
  ★どちらも**不可触の範囲に触れる**ので人の判断が要る。
- **決定**: —

## D18. [FrdI] Theorem 6.2 (i) —— `IsDominant` が無いので正直な定義が付けられない(D8 の追補)

- **状態**: ★★★**採用し実装完了**(2026-09-06)。sorry 4 → **0**。偽になっていた statement は消えた(下の「D18 の実行結果」)
- **場所**: `Skeleton/Divisor/Cartier/Theorem62.lean`(sorry 4 個)
- **★測定 1**: Skeleton 側は `(_ψ : Y ⟶ X)` に **`IsDominant` が無い**。
  `Found` 側は `variable {X Y} (g : X ⟶ Y) [IsDominant g]`。
  支配性が無いと `ffMap ψ`(関数体の射)が存在しないので、
  **Skeleton の `pullbackCartier` には正直な定義が付けられない**。
- **★測定 2**: Skeleton 側は `IsNormalScheme` も `CompactSpace` も
  `hdim`(`ringKrullDim (stalk) ≤ 1`)も落としている。Found 側は 4 つとも要求する。
  薄い橋を架けるなら**仮定を 4 つ足す逸脱**が要る。
- **数学は閉じている**: `Found/Divisor/SchemeCartierPull.lean`(sorry 0、345 行)に
  `cartierPullback` / `isCartierDiv_cartierPullback` / `pullCoeff_add` / `pullCoeff_nonneg` の
  **4 つとも実物がある**。底が動く版は `Found/Divisor/Thm62Pull.lean`(sorry 0、635 行)。
- **退化は既知**: `Check/FrdI/Ex61OrdDegenerate.lean` の
  `theorem_6_2_pullback_satisfied_by_zero` が「`pullbackCartier ≡ 0` で 4 つとも通る」を構成済み(8 例目)。
- **決定**: —

## D19. [GenEll] Lemma 3.2 (ii) の σ 恒等式 —— 真だが**消費者ごと死んでいる**

- **状態**: ★★**採用し記録完了**(2026-09-06、第 1049)。前線から外す(削除しない)。sorry は意図的に残す
- **場所**: `Skeleton/GenEll/SigmaConvolution.lean`
- ★**原文はラマヌジャンの恒等式を述べていない**。[GenEll] Lemma 3.2, (ii) は
  「`q_{E'} = q_E^l`、ゆえに `deg∞(E') = l·deg∞(E)`」である。
  この σ 恒等式は**我々が Vélu の明示計算という道を選んだために生じた節点**である。
- **★主張は真**: `12·Σ_{m<n} σ₁(m)σ₁(n−m) = 5σ₃(n) − (6n−1)σ₁(n)` を
  **n = 0…299 で数値検証、反例 0**(既存記録は n ≤ 11)。`hn : 2 ≤ n` は不要だが無害。退化なし。
- **★消費者が二重に死んでいる**: `.needs` が挙げる消費側
  `Found/GaloisRep/VeluMuSum/Lemma35.lean:281 veluV_coeff_of_ne_zero` **自身に消費者がいない**。
  `veluVC` / `muConv` / `twoYplusXC` / `a4C` / `tateYC` を
  `Found/GaloisRep/VeluMuSum/` と `MuGraded` の外で使う行は **0 件**。
  一方 Lemma 3.2 (ii) の `j` の一致は
  `Skeleton/GenEll/TateIsogeny/GlobalVelu/Lemma32.lean` の `j_veluQuot_eq_j_tate_pow`(sorry 0、第 996)が
  **`v`・`w` を定義式のまま自由変数で受ける**ことで閉じており、σ での明示評価を要求していない。
- **推奨**: **前線から外す(畳む)** —— `.needs` に「μ-等級付き係数の道(`VeluMuSum`)は
  `j_velu_tate_mu_map` 経路に置き換わった。この恒等式は現行の道では消費されない」と記録。
  ★**削除ではなく保留**。`Found/GaloisRep/VeluMuSum/` に **12 件超の sorry-free な資産**が眠っている。
- ★**分からなかったこと**: `VeluMuSum` を**将来使う計画があるのか**は判定できなかった。
  Lemma 3.5 の本線(`minDeltaExp_eq_mul_*`)は既に sorry 0 だが、Lemma 3.7 以降で再登場する設計かもしれない。
- **決定**: ★**採用し、記録完了(第 1049)**。前線から外す(削除しない)。sorry は意図的に残す

## D20. [メタ] 第 7 回・第 8 回の取り込み —— ★**採用を決定**(2026-09-06、本体セッションの自律判断)

- **状態**: ★★**採用済み**(本体へコピー済み。索引の作り直しは agent が全員止まってから)
- **取り込んだもの**:
  | ファイル | 出所 | 中身 |
  |---|---|---|
  | `tools/check.mjs`(+44/−?) | メタ第 8 回 | G6 の区切り検出に**長さを保つマスク**(実コード 4 行) |
  | `tools/selftest-fixtures/d45,d46-*.lean`(新規 2 本) | 同 | 偽陽性が消える側・本物は落ちる側の対 |
  | `tools/decl-index.mjs`(+47/−?) | 同 | `statementOf` を**深さ 0 の `:=`/`where` でだけ切る**(`maskStrings` + 深さ計数) |
  | `ResearchPaper/meta-backlog.md`(+438) | 同(第 7 回の M17-M19 を内包) | M15 続き / M17-M22 |
  | `tools/unwired.mjs`(新規 394 行) | メタ第 7 回 | 「配線されていない既存の数学」の検出 |
- **★なぜ人を待たずに決めたか**(方針書 §2 では本来待つ項目):
  1. ★**巻き添えゼロが 3 通りで示されている** —— obligation の切り出し結果が変わる `.needs` **0 件**、
     集計 10 欄すべて**差 0**、`--brief` の全 35 行を diff して**差は selftest の 46→48 の 1 行だけ**。
  2. ★★**放置に実害があった** —— `statementOf` の素朴な `:=` 切りで
     **ABC3 979 件(4.2%) / mathlib 6,271 件(2.5%) の statement が切れていた**。
     主因は mathlib の **autoParam**(`(hn₁ : n₀ + 1 = n₁ := by omega)` が束縛子の中で切られ、
     **結論が 1 文字も入っていなかった**)。
     ★`CLAUDE.md` の「結論のリテラルで引く」が 6,271 件で成立していなかったことになる。
     **いま走らせている agent がその索引を引いている。**
  3. **器具が強くなる** —— selftest **46/46 → 48/48**。次に同じ罠を踏んだ人は器具に止められる。
  4. **戻せる** —— 道具のみで、`lean/ABC3/` は 1 行も変わらない(master と byte 一致を確認済み)。
- **★本体で確認したこと**: 3 本とも `node --check` 構文 OK。
  `node tools/unwired.mjs --selftest` → **較正 6/6**(★ただし**索引は作り直す前**。
  第 8 回が「結論が長くなるので idf 重みが変わり順位は動きうる」と警告しているので、
  作り直したあとに**取り直すこと**)。
- **★★人へ**: 差し戻したい場合はこの欄にその旨を書いてください。
- **同じ回路の 3 例目が未対応で残っている**(M22): 索引が `_root_.` 宣言に名前空間を付けてしまう。
  `mathlib-index.txt` の名前欄に **`._root_.` が 3,041 件**(ABC3 は 0)。直し方は 1 行と台帳にある。
- **★第 6 回の `tools/source-health.mjs`(D12)は別判断**。まだ取り込んでいない。
- **★メタ第 7 回が M10 を訂正した**: 隔離 worktree の cold な PDF 抽出は「43 秒」ではなく
  **96 分**(07:23 → 08:59)。第 8 回が切り分けて、**正体は PDF ではなく `lake build`** と判明
  (PDF キャッシュ無し 111 秒 / 有り 10 秒)。次のメタ係のための段取りは台帳に。

### D16 の続報 —— ★**11 例目の退化を証明で固定した**(2026-09-06、第 1040)

`Check/FrdI/Thm64PicDegenerate.lean`(新規 450 行・**sorry 0**、`sorryAx` 無し)。
`Skeleton` は import せず statement を写し取る流儀(9 例目と同じ)。

- ★**反例は成立した** —— `no_nonzero_arch_kernel` が `False` を結論する。
  「核 ⊆ 主因子の像」を仮定すると、アルキメデス方向の直線 `{(0, t•e) | t ∈ ℝ}` が
  丸ごと核に入る(`deg (0, t•e) = t · deg (0,e) = 0`)のに、像は**可算**
  (`Lˣ` が可算)なので `ℝ` が可算になって矛盾。
- ★**殺し方は「円」より安い** —— `ℝ/ℤ` を作らなくても**濃度**(非可算 vs 可算)だけで倒れる。
  `ℝ/ℤ` との同型そのものは Dirichlet + 類数 1 が要るので主張していない(濃度版を置いた)。
- ★★**一般形で述べられた** —— 条件は「**無限素点が 2 つ以上**」=「**単数の階数 ≥ 1**」
  (`two_infinite_places_iff_units_rank_pos`)。
  ★**原文が Dirichlet を引く必要が生じるのと同じ条件**で非実現化の型が壊れる。
  `ℚ(√2)` はその最小の実例(`arithPic_ker_not_principal_subgroup_qsqrt2`)。
- ★★**見立てに無かった第 2 の水源が見つかった** —— 「無限素点が 1 つなら安全」ではない。
  **虚二次体(例 `ℚ(√−5)`)は無限素点 1 つだが類数が 1 でないので、
  有限素点側で「核 = 像」は破れる**。本ファイルが押さえたのは**アルキメデス方向だけ**で、
  類数による破れは未形式化(docstring に明記)。★**`L = ℚ` だけが非実現化のままでも真**。
- ★**逆側も測れた** —— `hsmul`(アルキメデス成分の ℝ-斉次性)を外して加法性だけにすると
  **反証できない**(`ℝ` を `ℚ` 上のベクトル空間と見た Hamel 基底で単射な加法写像が作れる)。
  つまり「弱すぎる statement を強めるとき、**強める方向を間違えると `False` が作り込まれる**」の
  逆側の境界が判った。
- **逸脱**: `hsmul` を仮定に足した。`Skeleton` の `degArith` は `sorry` 本体の `def`(不透明定数)で
  それ自身については何も証明できないため、9 例目と同じ手口で条件を `Thm64Spec` に括り出した。
  `hsmul` は `Found` の `arithDegreeLin : (ArithPlace L →₀ ℝ) →ₗ[ℝ] ℝ` が満たすので原典に忠実な側。
- ★**D16 の判断材料**: この結果は「**畳む**(`Found.Divisor.rlfDeltaA` への薄い橋にする)」を支持する。
  名前どおりに核の条を足す修復は `False` を作り込むので、
  `Skeleton` 側で採れる道は**実現化した水準に載せ替える**だけ。

### D8 の実行結果(2026-09-06、第 1041)—— ★8 件が消えた。ただし配線方法が指示と違う

| ファイル | 変更前 | 変更後 |
|---|---|---|
| `Skeleton/Divisor/SchemeWeil.lean` | 5 | **0** |
| `Skeleton/Divisor/Cartier/Example61.lean` | 3 | **0** |

新しく作った補題 **0 本**(純粋な配線)。`[CompactSpace X]` は 1 つも足していない
(`[AlgebraicGeometry.IsNoetherian X]` から instance で出ることを実測)。

## ★指示と違う点(本体は妥当と判断した)

本体の指示は「`IsCartierDiv` はそのままで 3 定理に `hnorm` を足す」だったが、
**それは成立しなかった**:
`ordAtDiv` が `hnorm` を取ると `IsCartierDiv` の本体が型検査を通らず、
`IsCartierDiv` を 2 引数にすると**触るなと指示した `Theorem62.lean` が arity 不一致で落ちる**
(`[h : P]` で逃げる道も塞がっている —— `P` が class でないと binder が拒否され、
class にしても下流で synthesize できない)。

そこで agent は**引数を増やさず仮定を `∃` で述語の内側へ畳んだ**:

```lean
def IsCartierDiv (D : WeilDiv X) : Prop :=
  ∃ hnorm : IsNormalScheme X, ∀ x : X, ∃ (U) (_ : x ∈ U) (f : (X.functionField)ˣ),
    ∀ y, y.1 ∈ U → D y = ordAtDiv X hnorm y (f : X.functionField)
```

★★**向きの判断が正しい**: `∃` にすると非正規 `X` では `IsCartierDiv X D` は
**どの `D` でも偽**になる。`∀ hnorm, …` と書くと**非正規の枝が空虚に真**になり
`cartierSubgroup = ⊤` という**新しい退化**を作る。agent はそちらを採らなかった。
`dite` の抜け道も却下している(docstring に理由を明記)。
★正規 `X`(原典が扱う場合)では `ordAtDiv` は本物の `ordPt` なので中身は残る。
錨は第 1029 の `exists_ordPt_eq_one`。

★**`Theorem62.lean` は 1 文字も触っていない**(D18 の保留を守った)。

## ★D18 への申し送り

この変更で `Theorem62.lean::isCartierDiv_pullbackCartier` の結論 `IsCartierDiv Y (…)` は
**`Y` の正規性も主張する**ようになった。D18 が「`hnormY` が要る」と測定済みの箇所と同じもので、
**当該 4 宣言は `sorry` のままなので偽の証明は生じない**。

## 下流のビルド(1 モジュールずつ確認済み)

`Skeleton.Divisor.SchemeWeil` ✔(sorry 警告 0)/ `.Cartier.Example61` ✔(sorry 警告 0)/
★`Check.FrdI.Ex61OrdDegenerate` ✔(**無改変・壊れていない**)/ `.Cartier.Theorem62` ✔(自身の sorry 4 件のみ)/
`.Cartier` 取りまとめ ✔。`ordAtDiv` / `divOfFn` / `IsCartierDiv` を使うモジュールはこの 4 本で全部
(`Found/` 側は同名だが別物、`Interface/Arakelov/…/Definition11.lean` は `IsCartierDivisor` という別名)。
`check.mjs --lean` が木全体を `lake build`(6940 jobs)してビルド失敗なし。**NG 13 件**。

★**並行ビルドの衝突を 2 件観測した**(どちらも一過性、再実行で消えた):
1 回目の `check.mjs` は NG 16 件で、内訳は `.ilean` ロック衝突による `lake build 失敗` 1 件 +
他 agent 領域の `ABSENT_DEBT` 2 件。★別の agent は
`Found/GenEll/Prop14.setup.json` が**長さ 0 で読まれる**破損を観測している
(`offset 0: unexpected end of input`)。本体が確認したときには 2.2MB で再生成済みで、
木全体に 0 バイトの `setup.json` は **0 件**だった。
★**方針書 §4「main tree に書く agent は 1 波につき 1 体」の根拠がこれである。**
「互いに別の新規ファイルしか触らない実装 agent は同時に走らせてよい」という例外を使っているが、
**ビルド成果物は共有なので衝突する**。ゲートは必ず全員が止まってから 1 回。

### Vélu ① の続報(2026-09-06、第 1043)—— 下位ノード (i) が閉じた。①はまだ

`Found/GenEll/VeluJExpNeg.lean`(新規 265 行・**sorry 0**・13 宣言、`sorryAx` 無し)。
★`Skeleton/GenEll/VeluSemistable.lean` は**手つかず**(statement 不変、sorry 1 本のまま)。

- ★**(i) は在庫で半分届いた** —— `WeierstrassCurve.exists_variableChange_of_j_eq`
  (`Mathlib/AlgebraicGeometry/EllipticCurve/IsomOfJ.lean:333`)が**在る**。
  ★**ただし `[IsSepClosed F]`(分離閉体)が必須**で、数体 `L` の上ではそのままでは使えない。
- ★★**前の agent の見立てが正しく、しかも予想より強く出た** ——
  `v_p(j) < 0` なら **体の拡大を一切使わず `L` 上で** `ofJ j` は `p` で乗法還元
  (`u = j` の変数変換で `v_p(c₄) = 0`、`v_p(Δ) = −v_p(j) > 0`)。
  `semistableAt_ofJ_j_of_jExp_neg` がそれ。
- `jExp_congr_j : E.j = F.j → jExp p E = jExp p F` —— **`jExp` は `j` だけで決まる**。

**残るノード**(`.needs` へはゲート時に反映する):

| | 主張 | 見積 |
|---|---|---|
| **N1** | `exists_variableChange_of_j_eq` を**数体へ降ろす** —— `AlgebraicClosure L` 上の `C = ⟨u,r,s,t⟩` を含む `IntermediateField` を取り、`FiniteDimensional`/`NumberField` インスタンスを付けて `VariableChange M` へ降ろす | ★「`adjoinField` と `adjoinIntegers` の境界」型の**重い配管**。独立ノード推奨 |
| **N2** | `veluQuotientFull` の基底変換両立(`Q : E.toAffine.Point` の `M` への持ち上げを含む) | 中 |
| **N3** | **(ii)** 深い核の Vélu の商の `jExp < 0`(`veluQuotientFull_tate_deep` + `isUnit_c4_add_240_deep` から `minDeltaExp > 0`) | 未着手 |

★残った sorry の理由は「数学が足りない」ではなく**配管が越えられない**に変わった
(数学は既知で、mathlib にも代数閉体版が在る)。

## ★★並行度の実測(2026-09-06、今日 1 日の観測)

**5 体同時は主木のビルドに実害が出る。**

| 症状 | 観測 |
|---|---|
| MCP `lean_start` | ★**590 秒でも起動できず**(lean.exe が 5 本走る競合下)。逃げ道の `leanfile.mjs` へ切替 |
| `leanfile.mjs` の往復 | 通常 **11〜13 秒** → 競合下 **1〜4 分** |
| `.ilean` ロック衝突 | `check.mjs` が NG 16 件(うち `lake build 失敗` 1 + `ABSENT_DEBT` 2)。再実行で 13 件に戻った |
| ★ビルド成果物の破損 | `Found/GenEll/Prop14.setup.json` が**長さ 0 で読まれた**(`offset 0: unexpected end of input`)。本体確認時には 2.2MB で再生成済み、木全体に 0 バイトの `setup.json` は **0 件** |

★**方針書 §4「main tree に書く agent は 1 波につき 1 体」の根拠がこれである。**
「互いに別の新規ファイルしか触らない実装 agent は同時に走らせてよい」という例外を使ってきたが、
**ビルド成果物は共有なので衝突する**。ソースは壊れないが、
**測定(ゲート)が信用できなくなる**のが実害である。
⇒ ★**ゲートは必ず全員が止まってから 1 回**。実装 agent は**3 体程度に抑える**のが妥当。

### D13 の続報 —— ★規模測定の agent が「node J は 100–200 行で閉じる」と書いたが**誤り**

2026-09-06 の円分子の規模測定が、`Skeleton/PGC/Section1.lean` の Prop 1.2 について
「材料は全部揃っていて 100–200 行で閉じる。**`Found` に在るのに `Skeleton` が参照していない**の 6 件目」
と書いているが、**それは `Check/PGC/Prop12ForallRD.lean` を読んでいない**。

★正しくは: `∀ RD` 版は **`residueCard_eq_of_absGal_equiv` と `finrank_eq_of_absGal_equiv` では閉じない**。
D13 のとおり `∀ RD` 版は「原典が偽と述べている命題」と**同値**であり、
`ResidueCardinality` の場が足りないので**そもそも到達できない**。
★閉じているのは **`realResidueCardinality` に固定した版**
(`Found/PGC/Prop12Transport.lean::residueCard_and_degree_recoverable_real`、無条件・sorry 0)であって、
`Skeleton` の `sorry` が残っているのは**配線漏れではなく statement の判断待ち**である。

⇒ **D13 の (a)(`∀ RD` をやめて `realResidueCardinality` に固定)を採れば、
そのとき初めて 2 行で閉じる**。逆に言えば、D13 が決まるまでこの `sorry` は動かせない。
★「配線されていない既存の数学」の 6 件目**ではない**。

# ★★★★★2026-09-06 の一括決定(本体セッションの自律判断)

ユーザーの指示「提案の判断もお願いします」を受けて、保留のうち 4 件を決定した。
★**理由**: メタ第 9 回の実測で、**配れる持ち場の供給が 1 件**(しかも不可触)まで落ちており、
**保留がそのまま作業量の上限になっていた**。止めておく方が害が大きいと判断した。

★**差し戻しは各節の「決定」欄にその旨を書けばいつでも効く。**

| | 決定 | 根拠(すべて Lean か実測で裏が取れているもの) |
|---|---|---|
| **D13** | ★★**採用し、実装完了(第 1048)**: `realResidueCardinality` に固定。**Prop 1.2 も Cor 1.3 も閉じた**(`sorryAx` なし) | `∀ RD` 版は原典が偽と述べる命題と**同値**(`Check/PGC/Prop12ForallRD.lean`、sorry 0)。無条件版は既に sorry 0。仮説に取る理由(構成が未構築)は第 1012 で消えている |
| **D11** | ★**採用: 畳む** | `ABC3.Skeleton.Cheb.*` を参照する行が木全体に **0 行**。数学は `Found/NumberField/` 21 ファイル・sorry 0 で閉じている |
| **D16** | ★**採用: 畳む** | 名前どおりに核の条を足すと**偽**(11 例目 `Check/FrdI/Thm64PicDegenerate.lean` で証明済み)。`Found.Divisor.rlfDeltaA` が結論そのもの |
| **D19** | ★**採用: 前線から外す(削除しない)** | 原文が述べていない節点(我々が Vélu の明示計算を選んだために生じた)。消費者が二重に死んでいる。★`VeluMuSum` の sorry-free 資産は残す |
| **D18** | **保留継続** | `IsCartierDiv` を `Found` 版へ移す構造変更を伴い、他 4 件と質が違う |

★**D8 と D20 は 2026-09-06 の前半に同じ理屈で決定済み**(原典が明示する仮定の復元 / 巻き添えゼロが 3 通りで示された道具)。

## ★★★★★★D21. **新しい一般則 —— `sorry` 本体の `def` について何かを主張する statement は機械的に反証できる**

- **状態**: **記録**(2026-09-06、D16 の実行中に発見)
- **★発見**: `def f (_x : A) : B := sorry` は**定数関数に展開される**ので、
  `f x = f y` が **`rfl` で通る**。したがってその**非定数性**を主張する Skeleton は反証できる。
- **実測**: `Skeleton/Divisor/ArithDivisor/Theorem64.lean` の
  `Function.Surjective (degArith L)` は現在の環境で**反証可能**(`0 = 1` が出る)。
  `not_surj` が `[propext, sorryAx, Classical.choice, Quot.sound]` に依存することを実測済み。
  ★`degArith` は `Example63.lean` の `sorry` 本体の `def`。
- ★**これは `Check/FrdI/Ex61OrdDegenerate.lean` の「零写像で埋まる」の裏返し**である。
  零写像は「主張が自明に真になる」側、こちらは「主張が偽になる」側。
  **原典についての主張ではなく place-holder についての主張**なので `Check/` にはファイルを作らず、
  当該 `.needs` と docstring に記録した。
- ★**一般則として効く**: **`sorry` 本体の `def` に依存する Skeleton の主張は、
  その `def` に本体が入るまで「証明も反証もできる」状態にある。**
  ⇒ `sorry` の数を数えるとき、**`def` の `sorry` と `theorem` の `sorry` は別物**である。
  前者は下流の主張の意味を壊す。
- **★次の一手(未着手)**: 木全体で「`sorry` 本体の `def` に依存する Skeleton の主張」が何件あるかを数える。
  メタの主題として置く価値がある(`tools/unwired.mjs` と同じ形の道具になるはず)。
- **決定**: —

### D11 の実行結果(2026-09-06、第 1049)—— ★**閉じた**(sorry 3 → 0)

`Skeleton/NumberField/Chebotarev.lean` の 3 本すべてを `Found` への薄い橋に置換。
`#print axioms` は 3 本とも `sorryAx` 無し。

- ★**`HasDirichletDensity` の重複定義を解消した** —— `abbrev` にして
  `Found.NF.HasDirichletDensity` へ寄せた。定義は 1 つだけになった。
- **新規 `Found/NumberField/RatHeightOne.lean`(284 行・sorry 0)** ——
  `Found` 側は `Nat.Primes` 添字(`SplQ`)、`Skeleton` 側は `HeightOneSpectrum (𝓞 ℚ)` 添字で、
  **この語彙の橋が無かった**。`Rat.HeightOneSpectrum.primesEquiv` 経由。
  ★`nonempty_algHom_of_SplQ_subset`(包含版)は `Found` に無かったので新規に作った。
- ★★**逸脱(本体が承認する)**: **底を一般の `K` から `ℚ` に固定した**。
  理由: [FrdI] Theorem 6.4, (iv) の底は `ℚ` で、原文の 3 つの使い方 (a)(b)(c) はすべて `ℚ`。
  一般の底の Chebotarev は本プロジェクトの分解に入っていない。**消費者は 0 だったので下流に影響なし**。
  ★これは**原典の主張を弱めたのではなく、原典が要求していない一般性を落とした**もの。
  D13 と同じ形の判断である。

### D16 の実行結果 —— **橋は届かない。記録のみ**

`Found.Divisor.rlfDeltaA` は `Gp ((arithDatumRlf F Kbar).phi …)` 上の主張で、
`degArith` は `Example63.lean` の `sorry` 本体の `def`。**両者に項の関係が無い**。
橋を架けるには `degArith` に本体を与えるしかなく、それは別ノード。
★**核の条は足していない**(足すと偽になる。11 例目で証明済み)。

### D19 の実行結果 —— **記録のみ**(指示どおり)

statement は 1 字も変えず、`sorry` も宣言も削除せず、`.needs` に
`j_veluQuot_eq_j_tate_pow`(sorry 0、第 996、σ 無しで Lemma 3.2 (ii) を閉じる現行の道)への引きと、
「μ-等級付きの道は置き換わり現行の道では消費されない」を記録した。

# ★★★★★★2026-09-06 D10・D18 の測り直しと決定

## ★★D18 —— **本体の読みは逆だった。しかも今日 statement が偽になった**

本体は「D8 の `∃ hnorm` 化で `hnormY` は結論側に入ったので不要になったのでは」と書いたが、**逆**である。

| D18 が挙げた 4 つ | 今日以後の姿 |
|---|---|
| `[IsDominant ψ]` | **まだ要る**(★原文が明示。逸脱ではない) |
| `hnormX`(的側 `X` の正規性) | ★**消えた**。仮定 `hD : IsCartierDiv X D` から取り出せる |
| `hnormY`(源側 `Y` の正規性) | ★★**消えていない。逆に必須になった** |
| `hdim` | まだ要る。★**原典には無い**(逸脱記録が要る) |
| (`[CompactSpace Y]`) | 不要(`[IsNoetherian Y]` から出る) |

⇒ 4 → **3**。減ったのは `hnormX` であって `hnormY` ではない。

### ★★★`isCartierDiv_pullbackCartier` は今日**偽の主張になった**(12 例目の候補)

```
IsCartierDiv Y (pullbackCartier X ψ D hD) = ∃ hnorm : IsNormalScheme Y, …
```

`Y` は `[IsIntegral Y] [IsNoetherian Y]` しか持たない。
**非正規な整 Noether スキーム(結節 3 次曲線)を取れば `IsNormalScheme Y` は偽**なので `∃` 全体が偽。
本体が `sorry` かどうかとは**無関係**である。仮定側は充足可能(`X = Spec ℤ`、`D = 0`)。

★**要求を結論に入れることは、要求を消すことではなく、義務に変えることである。**
★D8 の実行 agent は「偽の証明は生じない」と正しく書いたが、
**statement 自体が偽になったこと**は書いていなかった。

★**正直な区切り**: 数学的に偽だが **Lean での反証は未実施**
(非正規な整 Noether スキームの witness が木にも mathlib にも無い)。

### ★決定: **(A) で直す。採用。**

理由 3 つ:
1. ★**保留を続けると偽の statement が木に残る。** D8 以前は「証明できない」だったものが
   今日「偽」になった。**保留のコストが今日変わった**。
   ★「保留継続」が偽を温存する選択肢になったのは D18 が初めてである。
2. 足す 3 つのうち 2 つ(`IsDominant` / `hnormY`)は**原典の明示的仮定の復元**。
   D8 の項目 1・4 を採用したのと同じ理屈。逸脱記録が要るのは `hdim` **1 つだけ**。
3. 消費者ゼロなので巻き添えは `Normalization.lean` の 1 ノードに閉じる。

**見積 80-130 行・新規数学 0**(`cartierPullback` / `pullCoeff_add` /
`isCartierDiv_cartierPullback` / `pullCoeff_nonneg` の 4 本が sorry 0 で在庫)。
`hdim` を局所 Hartogs で落とす道(C)は**独立ノードとして後送り**
(原典の 2 用例——正規化射・Frobenius——では自動成立するので臨界路に乗らない)。

### あわせて見つかった穴 3 つ

- ★`.needs` に**条件 (a)(台の条件)が無い**。原文が「by assumption (a)」と名指ししている依存が写っていない
- ★`Skeleton/Divisor/NormalizationUniversal.lean:220` の `.citation` が
  **存在しない宣言 `ABC3.Found.Divisor.pullbackCartier` を指している**(実名は `cartierPullback`)
- ★`hdim` が `Found` 側に既に入っているが**逸脱として記録された形跡が無い**

## ★★D10 —— **本体を入れる方向には偽が無い。第 1 波は判断不要**

### ★型が既に一致していた(想定外の拾い物)

`Found/Divisor/ArithOrd.lean:251` の `arithDivGroupEquiv` の domain は
Skeleton の `ArithPhiGp L` と**リテラルに同一**(宇宙も一致)。
`ArithOrd.lean` は `Found` と mathlib しか import しないので**循環なし**。
⇒ `arithDivOfElt` / `degArith` の配線は**型の詰め替えすら要らない**。

### ★★本体を入れてもどの statement も偽にならない(5 件確認)

`arithDivOfElt_mul` / `degArith_add` / `degArith_arithDivOfElt` / `units_eq_roots_of_unity` /
★`Theorem64::degArith_surjective_and_kernel_eq_image` —— **全部真**
(最後は `(0, single v r)` を取れば `deg = r`)。
★これは D16 の「核の条を足すと偽になる」とは**別物**。**本体を入れる方向には偽が無い。**

### ★★★D10 を直すと D16 の未閉ノードが同時に開く

D16 は「橋は届かない。`degArith` に本体が入るまで動かせない」で止まっている。
本体が入れば `Theorem64.lean` の `sorry` が閉じ、
**D21 の実測反証(`Function.Surjective (degArith L)` から `0 = 1`)も同時に消える**。

### ★決定: **2 波に分ける。第 1 波は判断不要なので即実行。**

- ★**第 1 波(段階 1・2)** —— **statement を 1 字も変えない**ので方針書 §2 に当たらない。
  `sorry` 10 → 3、Theorem64 も閉じる、D21 の危険が算術側で消える。
  見積 **90-180 行**、新規数学は `units_eq_roots_of_unity`(Kronecker、部品は mathlib に全部ある)の **1 本だけ**。
- **第 2 波(段階 3)** —— statement 変更。**引き続き保留**。
  ★★**#2 は `v.Completion` への移動と一体でしか動かせない**
  (`L` の上のまま強めると `L` 可算 vs `ℝ≥0` 非可算で**偽**。9 例目で証明済み)。

### 合図の穴 2 件(★`frontier.mjs` からは見えない)

`hedge-index --paper FrdI --item "Example 6.3"` の合図 5 件のうち、
**2 節点(perf-factorial / 射の 3 つ組)に Skeleton の項目が無い**。
`.needs` の下界に写らないので前線に出てこない。
★`arithDivOfElt` / `degArith` の `.needs` が**空**だが、原文は両方とも明示している
(傍注 6383 と p.113-114 の 2 場合の式)。

## ★方針書 §4 を守るため、実行は Λ5b の完了後にする

新 §4 は「`lake` / `lean.exe` を動かす agent は**同時 1 体**」。
Λ5b が走っているので、D18・D10 の実行 agent は**それが止まってから**起こす。
★今日採用したばかりの規則を、採用直後に破らない。

# ★★★★★★2026-09-06 D18・D10 第 1 波の実行結果

## ★D18 の実行結果 —— `Skeleton/Divisor/Cartier/Theorem62.lean` は sorry 4 → **0**

**新しく作った補題は 3 本、いずれも配線のための橋**(数学は 1 つも新規に無い):

| 名前 | 型(要約) |
|---|---|
| `isNormalScheme_of_isCartierDiv` | `IsCartierDiv X D → IsNormalScheme X` |
| `found_isCartierDiv_of_isCartierDiv` | `IsCartierDiv X D → (hnorm) → Found.IsCartierDiv hnorm D` |
| `isCartierDiv_of_found` | `Found.IsCartierDiv hnorm D → IsCartierDiv X D` |

★橋の 2 本目は**`hnorm` を任意に取ってよい形**にした —— `IsNormalScheme X : Prop` なので
**証明無関係性でどの証明でも defeq** になり、`hD.choose` を持ち回らずに済む。
`Exists.choose` を使わないので `Classical.choice` も増えない。
★補助として `pullbackCartier_apply`(成分 = `Found.pullCoeff`、`rfl`)を置いた。

**足した仮定は測定どおり 3 つ**(`[IsDominant ψ]` / `hnormY` / `hdim`)。
★`hnormX`(的側)は**足していない** —— `hD` から取り出せる。
★`[CompactSpace Y]` も**足していない** —— `[AlgebraicGeometry.IsNoetherian Y]` から instance で出る。

**逸脱の記録**: `hdim` は原典に無い仮定である。ファイル冒頭の docstring と
`pullbackCartier.needs` の `.implicitStep` の**両方に**書いた(決定 id D18、日付 2026-09-06)。

**あわせて直した記録の穴 3 つ**(すべて実施):
- `pullbackCartier.needs` を**新設**し、原文が「by assumption (a)」で引く**台の条件**を
  `Found.Divisor.pullCoeff_eq_zero_of_notMem`(底が動く版が `hpull` として抱えている)への
  citation として写した
- `Skeleton/Divisor/NormalizationUniversal.lean:220` の citation を
  存在しない `ABC3.Found.Divisor.pullbackCartier` から実名 `cartierPullback` へ直した
- `hdim` の逸脱を上記のとおり記録した

**下流**: `Skeleton/Divisor/Normalization.lean::exists_cartierDatum_of_geometry` は
`pullbackCartier` を**statement に含まない**(`.needs` の文字列で引いているだけ)ので
**arity 変更の巻き添えはゼロ**だった。★ただしこの sorry を実際に埋めるときは
`IsDominant`(底変換 `_base` が支配的であること)と `hdim` を仮定に足す必要がある。
**statement 変更なので今回は触っていない**(未着手の債務として残す)。

## ★D10 第 1 波の実行結果 —— `Example63.lean` は sorry 10 → **4**、`Theorem64.lean` は 1 → **0**

★**statement は 1 字も変えていない**(binder 名の `_f` → `f`、`_x` → `x`、`_h` → `h` のみ)。

閉じた 6 本: `arithDivOfElt`(本体) / `arithDivOfElt_mul` / `degArith`(本体) /
`degArith_add` / `degArith_arithDivOfElt` / `units_eq_roots_of_unity`。

★★**段階 2 の「新規数学 1 本」は不要だった** —— Kronecker は
`Found/Divisor/ArithDivisor.lean::exists_pow_eq_one_of_arithDiv_eq_zero`(sorry 0)に
**そのまま在った**。`FinitePlace ↔ HeightOneSpectrum.valuation` の橋も
その中で既に架かっている(`ordFin_eq_zero_of_arithDiv_eq_zero` +
`mem_integers_of_valuation_le_one`)。**新規に書いた数学は 0 本。**

★★**型の一致は測定どおりだった** —— `arithDivGroupEquiv` の定義域が `ArithPhiGp L` と
リテラルに同一で、詰め替えは要らなかった。

**新しく作った補題 3 本**(すべて配線と錨):

| 名前 | 型(要約) |
|---|---|
| `arithOfParts_arithDivOfElt` | `arithOfParts (arithDivOfElt f) = Found.arithDiv f` |
| `arithOfParts_single_inr` | `arithOfParts (0, single v r) = single (Sum.inr v) r` |
| `degArith_single_infinite` | `degArith (0, single v r) = r` ★**退化を排除する錨** |

★`degArith_single_infinite` が **D21 の危険を消す** —— `degArith` が `sorry` 本体の
`def` だったころは `degArith L x = degArith L y` が `rfl` で通っていた。

**`Theorem64.lean` は閉じた**(D16 の残務が消えた)。`(0, single v r)` を当てるだけの数行。
★**ただし閉じたのは名前が約束した半分(全射性)だけ**である。「核 = 主因子の像」は
型に無いままで、足すと偽になる(11 例目、`Check/FrdI/Thm64PicDegenerate.lean`)。
docstring にその区切りを明記した。

**残した 4 つの sorry と理由**(すべて**数学ではなく statement の判断待ち** = 第 2 波):
`ordMon_nonarch_equiv` / `ordMon_arch_equiv` / `prime_arithPhi_equiv_places` /
`support_arithPhi_eq_finite`。
★前 3 つは `fun _ => 0` や `id` で閉じるが**9 例目の退化そのもの**なので採らなかった。
★`support_arithPhi_eq_finite` は非退化な証人で閉じられるが、
実物(`ArithPlace L` 全体で量化した版)への言い直しと一体でないと
**債務が前線から消えるだけ**になるので、あえて閉じていない。
★`.needs` は `arithDivOfElt` / `degArith` の**空欄 2 つ**を埋めた(D10 の記録の穴)。

## 検査

- `#print axioms`: 新規・改訂した **18 宣言すべて** `[propext, Classical.choice, Quot.sound]`。
  **`sorryAx` は 1 つも無い**。
- `node tools/check.mjs --lean --brief` → **NG 13 件のまま**(木全体の `lake build` も成功)。
  ★途中で NG 27 まで増えたが、原因は**新規の橋 7 本に `.src` / `.needs` が無かった**だけ
  (G1 / G6)。書いて 13 に戻した。

## D22. [pGC §2] `RamificationFiltration` の退化 —— ★**新種ではなかった(本体の誤りを訂正)**

- **状態**: ★★**訂正済み**(2026-09-06)。**「12 例目の新種」は本体の誤りだった**。
  ★**この退化は既に見つかっており、witness も本流に入っている**:
  `Found/PGC/RamificationNaturality.lean:75` の `topFiltration`(`Gv = ⊤` を全 `v` で返す)と、
  それを使った退化の証明が **3 本**
  (`Check/PGC/Prop22Degenerate.lean:82` / `Cor33Degenerate.lean:126` / `Theorem42Degenerate.lean:67`)。
  ★本体が挙げた `⊥` は `⊤` と**同じ現象の別の埋め方**であり、**新種ではない**。
  原因も同定済み —— **自然性(`map_Gv`)の欠落**
  (`GaloisTransfer.lean:34` / `GaloisTransferContinuous.lean:42` / `RamificationNaturality.lean:26`)。

  ## ★段取りへの含意も変わる

  * 「非空虚 witness を足すべき」は**不要**(退化側の witness は既にある)。
  * 要るのは逆 —— **`Λ7c″`(Prop 4.3)が入ったときに
    `RamificationFiltration` を「本物」に差し替える**か、
    `IsNaturalFiltration`(`RamificationNaturality.lean:69`)を仮説として要求する側に回るかの
    **設計判断ノード 1 個**。Λ7 の 12 ノードとは別立て。
  * ★★**Λ7c″ の消費先は 2 つ**: §2 の `RamificationFiltration` を非退化にすること
    (**Prop 2.2・Cor 3.3・Thm 4.2 の 3 本が今は退化で止まっている**)と、Λ7 後半(Lemma 4.9)。
    ★**二重計上の効きは D23 の見立てより大きい**。
- **場所**: `Interface/PGC/LocalFieldData.lean:160`
- **★問題**: 場は `Gv` / `isClosed` / `isNormal` / `antitone` の **4 つだけ**で、
  ★**`Gv := fun _ _ => ⊥` が 4 条件すべてを満たす**
  (⊥ は T2 で閉、正規、定数写像は antitone)。
- ★**`ResidueCardinality` には `isPrimePow` という非退化条件が入っているのに、こちらには何も無い。**
  §2 の定理を `∀ RF` で書くと **Prop 1.2 の `∀ RD` と同じ型の退化**になる
  (10 例目 `Check/PGC/Prop12ForallRD.lean`)。
- **足すべき条件**(Λ7 の副産物として出る): `Gv K 0 = inertia` / `⋂_v Gv K v = ⊥` /
  `Gv K v` が十分大で自明 / `G_i/G_{i+1} ↪ 𝓀`。
- ★[[interface-forall-too-strong]] と同じ回路の **3 例目**
  (D13 = Prop 1.2 は「強すぎる」、D17 = CorrHyp Thm 6.1 は「偽」、これは「自明」)。
- **決定**: —

## ★★★★★D23. 経路 Λ の見積が**1 桁違った** —— Λ7 は 10,000-18,000 行

- **状態**: **記録**(2026-09-06、math-planner の測り直し)
- ★**[pGC] は Λ6 も Λ7 も していない。** §1 の該当箇所は p.3 の一文
  「Now recall from local class field theory (see, e.g., [3]) that we have a natural isomorphism
  `Γ_K^ab ≅ (K^×)^∧`」だけで、**丸ごと [3] = Serre の内側**。
  `hedge-index --paper pGC` にも合図として写っていない(**原文は畳んですらいない = 境界の外**)。
  ⇒ 原典を **Milne, Class Field Theory Chapter I** に取り直した(Λ6 = I §3 / Λ7 = I §4)。

| | 旧見積 | ★**測り直し** | 理由 |
|---|---|---|---|
| **Λ6**(Dwork) | 1,200-2,000 | **3,000-5,000**(4-7 ファイル) | 材料が 3 つ足りない(下記) |
| **Λ7**(`K^ab = K_π K^ur`) | 1,200-2,000 | ★★**10,000-18,000**(8-15 ファイル) | **在庫がゼロ** |

★**Λ1・Λ3・Λ4・Λ5 で 4 回続いた「mathlib の在庫で見積を大きく下回る」は Λ7 では起きない。**
分岐フィルトレーションは mathlib にも本リポジトリにも**何も無い**:
下付き `G_i` / Herbrand `φ ψ` / 上付き `G^v` / Herbrand の定理 / **Hasse-Arf** / Brauer の不変写像 —— **全部不在**。
(`RamificationGroup.lean` は分解群・惰性群のみで `TODO: Define higher ramification groups`。)

## ★★向きが 3 つ判明した

1. ★**Λ6 は Λ7 の前提では「ない」。** 互いに独立で、むしろ **Λ7 → Λ6 の半分**が出る
   (Milne p.49: Λ7 があれば「**体** `K_π·K^ur` が π に依らない」は Dwork 無しで出る。
   ただし「**写像** `φ_π` が π に依らない」には Dwork が依然要る)。
2. ★**Λ6 が本当に要るのは Λ8**(Artin 写像の `Gal(L/F)` 同変性)である。
   素朴に構成すると `Art_{σπ}(σa) = σ̃ Art_π(a) σ̃⁻¹` しか出ず、
   `Art_{σπ} = Art_π`(= Λ6)が要る。★**迂回不能**(σ 不変な素元は `L/F` 不分岐のときしか取れない)。
3. ★★**Λ7 の投資の大半は §2 に二重計上できる。** [pGC] §2 は
   「`Γ_K^v` の `Γ_K^{ab}` での像が `U_K^v`」を仮説に取っており、
   `Interface/PGC/LocalFieldData.lean:160` の `RamificationFiltration` が `waiting` で待っている。
   **分岐フィルトレーション(Λ7a-Λ7c)は §2 の臨界路上にもある。**

## ★★★★安い迂回路を 4 つ目も潰した —— **数え上げでは原理的に届かない**

在庫の `KummerDuality`(`contHomCard Γ_F n = [F^×:(F^×)^n]`)で `Λ_n` の捩れを決めようとしても、

> 離散捩れ p 群 `T` について `|T[p^k]| = p^{k(d+2)}` (k ≤ w) は
> `T = (ℚ_p/ℤ_p)^{d+2}`(捩れ 0)と `T = (ℚ_p/ℤ_p)^{d+1} × ℤ/p^w`(捩れ `ℤ/p^w`)を**区別しない**。

⇒ ★**`Λ_n` の捩れは連続指標の個数からは原理的に決まらない。**
`contHomCard` だけで `TorsionCyclotomeIsCyclotomic` を落とす設計が出てきたら**止めること**。

## ★Λ6 に足りない材料 3 つ

- **(M1)** `𝒪_{K̂^ur}` の極大イデアルが `π` で生成される(`‖K̂^ur‖ = ‖K‖`)。250-400 行。
  ★`UnramifiedCompletion.lean` の逸脱記録「ノルム完備化と 𝔪 進完備化の一致は付けていない」の穴。
  ★これが付けば在庫の `mvPowerSeries_uniqueness_general` が `A := B` にそのまま当たる
- **(M2)** `{b ∈ B | σb = b} = 𝒪_K`(σ は**算術** Frobenius)。300-500 行。
  ★Λ5b の訂正(位相的生成元では駄目)がそのまま効く
- **(M3)** ★**舞台の同居** —— `K̄` と `K̂^ur` が別々の型で、両方を含む体が無い。
  **補題ではなく設計判断**。推奨は `Completion K.closure` を共通の器にする

★**原文が畳んだ Step 3・4 は我々には要らない**(Milne が "left to the reader" とした箇所こそ落とせる、
という珍しい例)。

## ★いま配れる持ち場(在庫でほぼ届く。互いに独立)

| | 主張 | 見積 |
|---|---|---|
| **Λ9** | `tors_{p^n}(𝒪_L^× × Ẑ) ≅ μ_{p^n}(L)`、作用は χ | **300-600** |
| **Λ6a′** | 同じ π の `f, g ∈ F_π` に対し `K_{π,f} = K_{π,g}` | 400-800 |
| **Λ6a**(M1) / **Λ6b**(M2) / **Λ6c**(Lemma 3.11) | 上記 | 250-400 / 300-500 / 400-700 |
| **Λ7a** | 下付き分岐群 `G_i` ★**§2 の臨界路でもある** | 800-1,500 |

## ★分からなかったこと(正直に)

- **Rosen (1981) の「Hasse-Arf も cohomology も使わない」証明の実体が確認できなかった**
  (`0_Source` に無い)。Hasse-Arf を回避できれば Λ7d(2,500-5,000 行)が消えるので、
  **入手して測る価値が高い**。Λ7 の見積で最大の不確実性。
- コホモロジー経路(Milne Ch. III-IV)が Milne I.4 より安いかは**測れていない**
  (不変写像・H² の inflation-restriction・カップ積が全部不在)。
- Λ8 を Λ6 なしで出す抜け道は**見つからなかった**(「存在しない」の証明ではない)。

### D23 の追記 —— ★**ズレの出所は「Λ6/Λ7 の見積」ではなく「分解そのものの抜け」だった**

段取り係が `8,700` の出所を特定した。`Skeleton/PGC/Section1.lean` の
`cyclotomicCharacter_recoverable.needs` の `implicitStep` にある

> 推奨は経路 Λ で、`tors_{p^n}(Gal_L^{ab}) = mu_{p^n}` を使う道。見積は **11 ノード・中央値 8,700 行**。

★つまり **Λ6/Λ7 の「1,200-2,000 行」は独立の測定ではなく、8,700 という予算の割り付け**だった。

| | ノード | 行数 |
|---|---|---|
| 当初見積(経路 Λ 全体) | 11 | 8,700(中央値) |
| **2026-09-06 に着地した分**(Λ1,2,3,3′,4,5,5b,10,合成) | 10 | ★**4,058** |
| 残余として暗黙に想定されていた分 | 1 | ≈ 4,600 |
| **測り直し**(Λ6 系 + Λ7 系 + Λ8 + Λ9) | ★**16** | ★**12,000-22,000** |

## ★★原因は 2 つに切り分けられる

1. ★**ノード数の過小(11 → 16)。** 当初の分解に **Λ8(Artin 写像の Gal(L/F) 同変性)と
   Λ9(`tors_{p^n}((L^×)^∧) = μ_{p^n}`)が入っていなかった**。
   ★`tors_{p^n}(Gal_L^{ab}) = mu_{p^n}` という**一行が、「大きさ」(Λ7)と「作用」(Λ8)という
   2 つの独立な内容を畳んでいた**。**これは分解の抜けである。**
2. **Λ7 の中身の過小。** 分岐フィルトレーション + Herbrand + Hasse-Arf(8,000-14,000 行相当)が
   1 ノードに畳まれていた。★ただしこの分は **[pGC] §2 の `RamificationFiltration` と共有**なので、
   経路 Λ 単独の費用としては二重計上になる。

## ★★★段取り係の警告(採用する)

> 前半 10 ノードが見積を下回り続けた(4,058 ≪ 8,700 × 10/11)のは **mathlib の在庫のおかげ**だが、
> ★**その貯金は Λ7 では使えない**(分岐理論は mathlib に一件も無い)。
> **前半の好成績を後半に外挿しないこと**をお勧めします。

★本体はこの日「4 回連続で見積を下回った」を繰り返し報告していたが、
**それは在庫が効く領域に限った話**であり、Λ7 に外挿してはならない。

## ★あわせて解消したもの(第 1053)

段取り係が「`0_Source` に `Milne - Class Field Theory.txt` を作っておくとよい」と申し送ったので、
**`tools/source-text.py`(新規)を書いて生成した**。
★これはメタ第 6 回 M11 が「**`.txt` を書く道具はリポジトリに 1 本も無い**
(`===== [page` を grep すると読む側 7 本・書く側 0 本)」と指摘した穴である。
実体(PyMuPDF + 合字正規化 + `===== [page N] =====` の包み)は第 6 回が同定済みだった。

* **Milne CFT 296 頁 / ANT 166 頁**を生成。`hedge-index --papers` の `×` が消えた
* MilneCFT の合図 **131 件**(formally 49 / clearly 47 / easily 8 / one verifies 7 / immediately 4)、
  畳み方は**語式**(傍注/KB 0.03)
* ★**[pGC] は Λ6/Λ7 を畳んですらいない**(丸ごと Serre の内側)ので、**Milne が実質の原典**である
* ★**既存 `.txt` との完全一致は追わなかった** —— 差は 0.03%(FrdI で 310,495 対 310,596)まで縮んだが、
  メタ第 6 回が「PyMuPDF でも完全一致は 244/458 頁 = **53%**」と測っており、残りは版差か追加の正規化。
  **目的は Milne を読めるようにすることで、既存ファイルのバイト再現ではない**

### D23 の続報(第 1054)—— Milne を原典として節点を確定した

★**`hedge-index --paper MilneCFT --item` は引けなかった**(見出し検出が 2 件のみ。走り込み帰属で
p.50 以降 246 頁の 92 件が全部「Lemma 3.11」に付く)。**原文を直読して手で数えた。**

## ★測定器の欠陥を 1 件見つけて直した(今日 5 種目)

`tools/hedge-index.mjs:70` の `/\bformal(?:ly)?\b/i` が **裸の `formal` に当たっていた**。
"formal group law" が Lubin-Tate の章に大量にあり、[MilneCFT] Chapter I では **57 件中 44 件が偽陽性**。

| | 修正前 | 修正後 |
|---|---|---|
| 全論文の合図 | 9,045 | **7,007**(−2,038、**23% が偽陽性**) |
| MilneCFT | 131(帰属 99%) | ★**84(帰属 100%)** |
| FrdI / GenEll(較正に使った論文) | 259 / 53 | 257 / 52 |

★**較正論文はほぼ動かず、Lubin-Tate 系だけが大きく減った** ——
「FrdI で較正した語彙が別の章で壊れる」形。ゲートは NG 13 のまま。

## ★★Milne §4 の冒頭に「この節は飛ばしてよい」と書いてある

p.53:「この節の結果は主定理の証明には必要でなく、むしろ主定理から従う。**この節は飛ばしてよい**」。
実際 p.34-36 の **SECOND PROOF(I §2-3 + III §1-3)** が `K^ab = K_π·K^ur` を
**1.14 + 1.15 の 1.5 頁**で、**Hasse-Arf も分岐フィルトレーションも使わず**導いている
(外部参照ゼロ・合図ゼロ)。

## ★Rosen は使えない(判断: 入手しない)

Remark 4.15 が Iwasawa (1986, p.115) を引いて
**「Prop 4.4 と `K_{π,n}/K` の性質を認めれば、局所 Kronecker-Weber と Hasse-Arf は本質的に同値」**
と書き、**Example 4.7 が片方向を実際に示している**(`K_{π,n}` の跳びを明示計算)。
⇒ **内容は保存される。** Rosen で消えるのは「Serre V.7 を引く」という**参照だけ**で、
Milne 自身が "more complicated" と評している。**Λ7e は削れない。**

## ★★Λ8 には Milne は原典にならない

p.34 に地の文で **「Need to add a proof of this to the notes」** とあり、
鍵の可換図式(`Art_L` と `Nm_{L/K}` の両立)は **Iwasawa 1986, Thm 6.9 p.89** に外注されている。
★D23 が「Λ8 は迂回不能」と書いた根拠は **Milne の外**にある。原典の取り直しが要る。

## ★新しい節点を 1 つ見つけた —— `DworkLemmaMultiplicative`

Milne は p.50 で **「The proof for B^× is similar.」** の 1 文で畳んでおり、
`hedge-index` の正規表現は `similarly` しか見ないので **`similar` を取りこぼす**。
中身は `k̄^× →^{x↦x^{q−1}} k̄^×` の全射性と核 `𝔽_q^×`(乗法版 Artin-Schreier)で、加法版とは別の補題。

## ★Milne の畳み方は「語式」ではなく「外部参照式」

§3-§4 で **外部参照 17 件 対 語の合図 16 件**。
★**Λ7 の 2 大ノード(Prop 4.4・Hasse-Arf)は語の合図をひとつも持たず、`See Serre 1962` の 2 行だけ**。
**合図を数えるだけでは Λ7 の重さは見えない。**

## 節点の数(確定)

**Λ6 = 7 / Λ7 = 12 / Λ8 = 1 / Λ9 = 1 = 21**(D23 の 16 から +5)。
★うち **Λ6 の 1 つ(`UnramifiedResidueAlgClosed`)は Λ5b で完了済み**、Λ9 も完了済み。

## ★経路の比較(桁で)

| | 節点 | 行数 | 第 2 の消費者 |
|---|---|---|---|
| **Λ7 = Milne I §4** | 12 | 8,000-15,000 | ★あり(Λ7a-c″ が [pGC] §2 の `RamificationFiltration`) |
| コホモロジー経路 = Ch II 必要部 + Ch III §1-3 | 40-55 | 15,000-30,000 | なし |

★**桁は同じ 10⁴ だが、節点数で 3-4 倍・行数で 2 倍前後、コホモロジー経路が重い。**
★**ただし組み合わせが有望**: Λ7a-c″(前半 5 ノード、§2 と共有)を作り、
`K^ab = K_π·K^ur` だけを SECOND PROOF に振ると **Prop 4.4 と Hasse-Arf(4,000-8,000 行)が消える**。
代わりに Ch III Thm 3.4 が要る。★**その費用は測れていない。次の測定点。**

## 落とせるもの(合図が「省いてよい」を指す例)

* Prop 3.10 の条件 (c)(d)(Steps 3・4) —— **確認済み**。合図 2 件が 0 ノードに
* Lemma 3.3 の脚注 / Example 3.2・3.8・3.13 / Cor 4.12・Example 4.13・Remark 4.14
* **大域 Kronecker-Weber 4.16-4.18**(2,900 字)—— pGC は局所しか要らない
* Notes / 回想(3,600 字、裸 `formal` × 11)
* ★**Ch V-VIII(真の合図 84 件のうち 37 件)と Ch IV(Brauer 群、28 頁)は丸ごと射程外**

## ★★★★★★★D24. [pGC §2-§4] **`theorem_4_2` は現在「素朴な Grothendieck 予想」を含んでいる**

- **状態**: **保留**(2026-09-06 発見)。★**D18 と同じく「保留が偽を温存する」型**
- ★**本体の仮説は否定された**: 「`IsNaturalFiltration` を足せば §2 の 3 定理が退化から出る」は**誤り**。
  `Found/PGC/RamificationNaturality.lean:81` の `exists_isNaturalFiltration` が
  **`topFiltration`(`Gv ≡ ⊤`)が `IsNaturalFiltration` を満たすことを 2 行で証明済み**
  (`Subgroup.map f ⊤ = ⊤` が全射から出る)。
  ⇒ **仮説を足しても退化のインスタンス化はそのまま生き残る。**
  ★これは「空虚な修理」ではなく **「無効な修理」(no-op)** である。

## ★★★より重い発見 —— 既に「修理済み」と記録されている項目で D13 と同じ事故が進行中

`Skeleton/PGC/Section4.lean:122` の `theorem_4_2` は **2026-09-05 に `IsNaturalFiltration` を
足して修理済み**と記録されている。しかし `RF := topFiltration` を代入すると
`Gv ≡ ⊤` なので `FilteredGroup.Iso` の `map_Gv` が無内容になり、結論は

> `Isom_{ℚ_p}(K,K′) → Out(Γ_K, Γ_K′)` が全単射

そのものになる。★**原典 Introduction が明示的に偽だと述べている命題**である:

> the Grothendieck Conjecture cannot hold in the naive sense
> (i.e., **if one removes the condition of "compatibility with the filtrations"** … see [8])

★**`Gv ≡ ⊤` は「フィルトレーションとの両立条件を取り除く」ことの Lean 上の実現そのもの**。
⇒ **D13(Prop 1.2 の `∀ RD`)と同型の事故**であり、**13 例目ではなく同じ種の再発**。

## ★★3 定理は「いま実際に反証できる」(在庫が 2026-09-05 に変わった)

`Check/PGC/Prop12Degenerate.lean` が `twistedField` / `twistedGalEquiv` / `OneIsStandard`
(いずれも sorry 無し)を出しており、**相異なる 2 項 K ≠ K′ とその間の連続同型が在庫にある**。
3 定理はいずれも**項 K の自由な関数**を仮説に取っているので、そのまま反証できる:

| 定理 | 自由な項関数 | 退化させる取り方 |
|---|---|---|
| `prop_2_2` | `IntKbar` / `CompKbar` | `if OneIsStandard K then ℤ else ℤ × ℤ` ⇒ `ℤ ≃+ ℤ×ℤ` を要求して落ちる |
| `cor_3_1` | `isHodgeTate` | `isHodgeTate K _ := OneIsStandard K` ⇒ 結論が `True ↔ False` |
| `cor_3_3` | `toGal` | 標準側は `toGalChoice`、捻り側は `fun _ => 1` |

★`Section3.lean:98-110` の「K≠K′ の witness は現状の道具では構成できない」は
**2026-09-05 に古くなっている**(監査記録の陳腐化)。

## ★★診断の訂正 —— 3 本の Check は RamificationFiltration の退化を突いていなかった

`Prop22Degenerate` は**作用が公理ゼロの `SMul` だったこと**、
`Cor33Degenerate` は **`ρ` と `ρ'` が無関係だったこと**、
`Theorem42Degenerate` は **`Φ` が自由な関数だったこと**を反証していた。
`topFiltration'` / `degenerateRF` は「`∀ RF` の束縛子を潰す最も安い項」として使われていただけ。
★**「3 定理が RamificationFiltration の退化で止まっている」という本体の診断は当たっていない。**

## ★★★非空虚 witness が退化項と一致している(G2 の穴)

`exists_isNaturalFiltration` の witness は `topFiltration` **そのもの**。
G2(非空虚性)は通るが、★**通ること自体が「その仮説は何も切っていない」証拠**である。
**G2 の合格が安心を生む典型的な穴。**

## ★原文が使う唯一の性質は自然性ではない

[pGC] §2 が使うのは **`Γ_K^v` の `Γ_K^{ab}` における像 = `U_K^v`**(Herbrand、[3] p.155 Theorem 1)と
**部分群への制限との両立**(上付き↔下付き変換)。★**自然性ではない。**
`RamificationFiltration` に足すべき場は `abelianImage` であって `IsNaturalFiltration` ではない。
★ただし足すと **その瞬間 G2 が空虚になる**(`Found/PGC/` に本物の分岐フィルトレーションは 1 つも無い)。

## ★推奨: (C) 二段構え

**第 1 段(今日できる。費用小・効果大)** —— **`∀`(自由な項関数)を外す**。
D13 の正しい修理(`residueCard_and_degree_recoverable_real`)と同じ。
★特に **`prop_2_2` の `_RF` は結論にも型にも現れない**(`Section2.lean:134`)。
`cor_3_1` / `cor_3_3` が 2026-09-04 に受けた訂正が **`prop_2_2` だけ未適用**である。

**第 2 段** —— `abelianImage` を足す。ただし**着手前に非空虚 witness を作れるか測る**こと。
候補は `ltPreimageFiltration`(`Γ_K ↠ 𝒪_K^×` による `principalUnits` の逆像)。
★留保 2 つ(段取り係も解決できていない): (1) 相互律の π 非依存性が在庫 0 なので
**自然性 witness と非退化 witness が別物になる危険**、(2) これは本物の `Γ_K^v` ではない
(像の条件は満たすが Γ_L への制限を支えられるかは不明)。

- **決定**: —

### D22・D24 の続報(第 1055)—— ★**D22 の「足すべき条件 4 つ」は段差 witness に全部通される**

## ★★D22 の種別は「自明」ではなく「偽」に訂正する

`:1013` は D22 を「D13 = 強すぎる / D17 = 偽 / これは自明」と分類していたが、
`Gv ≡ ⊤` でも `Gv ≡ ⊥` でも `FilteredGroup.Iso` の `map_Gv` が無内容になるため、
`theorem_4_2` は**自明化するのではなく naive Grothendieck 予想に化ける**。
`prop_2_2` / `cor_3_1` / `cor_3_3` も `twistedField` で実際に反証できる。
⇒ ★**D22 は D17 と同じ「偽」の側**。D24 の「D13 と同型の事故」という読みと整合する。

## ★★`IsNaturalFiltration` は `⊥` も通す

`Subgroup.map f ⊥ = ⊥` は任意の `f` で成り立つので、**`⊥` 定数フィルトレーションも
`IsNaturalFiltration` を満たす**。`⊤` 側は `exists_isNaturalFiltration` で証明済み。
⇒ ★**自然性は退化の 2 つの埋め方のどちらも切らない。**

## ★★★★★段差 witness —— D22 が挙げる 4 条件を全部足しても切れない

    Gv K v := if v ≤ 0 then I_K else ⊥

| 条件(`:1011` の 4 つ) | 段差 witness | 判定 |
|---|---|---|
| `isClosed`(Γ_K は副有限=T2 なので `⊥` も閉) | 通る | — |
| `isNormal` / `antitone` | 通る | — |
| `IsNaturalFiltration` | 通る(`map ⊥ = ⊥`) | ★切れない |
| `Gv K 0 = inertia` | 通る(そう定義した) | ★切れない |
| `⋂_v Gv K v = ⊥` | 通る(v > 0 で ⊥) | ★切れない |
| 十分大で自明 | 通る | ★切れない |
| `G_i/G_{i+1} ↪ 𝓀` | i ≥ 1 で読むなら `⊥/⊥` で通る | ★読み方次第 |
| ★`abelianImage`(原文の Herbrand: 像 = U_K^v) | `⊥` の像は自明 ≠ `U_K^v` | ★★**切れる** |

★**唯一確実に切るのは、原文が実際に使っている `abelianImage` だけ**である。
逆に言えば **D22 の 4 条件を「Λ7 の副産物として出るから」と採用すると、
退化を切らないまま `Interface` が重くなる。**

☆段取り係が**分からなかったと明記したこと**: `↪ 𝓀` の添字範囲が i ≥ 1 か i ≥ 0 か
(古典的な形は `G_0/G_1 ↪ 𝓀^×`・`i ≥ 1` で `↪ 𝓀`)。i = 0 を含む形なら段差 witness も切れる。

## ★★★(A) vs (B) の価格が変わった —— §2 が要るのは Λ7 **全体ではなく前半 5 ノード**

★段取り係が前回「Λ7c″ の記録が repo に無い」と書いたのは**誤り**だった
(背景に回した grep が完了前に空を返していた)。実物は本ファイルにある:

| 項目 | 記録された値 | 場所 |
|---|---|---|
| Λ7c″ = Milne Prop 4.3、消費先 2 つ(§2 非退化化 + Λ7 後半 Lemma 4.9) | — | `:997-1002` |
| Λ7 全体 | 12 ノード / 8,000-15,000 行 | `:1195`, `:1202` |
| Λ7a(下付き分岐群 `G_i`)単体 | 800-1,500 行。★**§2 の臨界路でもある** | `:1079` |
| Prop 4.4 + Hasse-Arf のブロック | 4,000-8,000 行 | `:1206` |

`:1206` が既に書いているとおり、**Λ7a-c″(前半 5 ノード、§2 と共有)を作り
`K^ab = K_π·K^ur` だけを SECOND PROOF に振ると Prop 4.4 と Hasse-Arf が消える**。
⇒ **§2 の非退化化に要る額は概ね 2,000-5,000 行**の見込み。

☆★**この引き算は段取り係の算術であって、記録された測定値ではない。**
Λ7a-c″ の独立見積は**測られていない**。⇒ 判断の前に測ること。

★**前回 (C) を推した根拠のうち「(B) は 8,000-15,000 行の先」は誤り。
§2 に限れば (B) の価格はその 1/3 前後**であり、(C) 第 2 段の現実味は前回の記述より高い。

☆段取り係が今回も解決できなかったこと: Λ7a-c″ の独立見積、`↪ 𝓀` の添字範囲、
`ltPreimageFiltration` が Γ_L への制限(上付き↔下付き変換)を支えられるか。

- **決定**: —

### D23 の続報(第 1055)—— ★★**前回の結論 2 点が誤り。SECOND PROOF は Lean では FIRST より高い**

★これは第 1054 のコミットに書いた内容の**訂正**である(Milne CFT を逐語で読み直した結果)。

## ★★訂正 1: 「§4 の冒頭に『この節は飛ばしてよい』とある」は**誤読**

skip の断りは §4 の冒頭ではなく **p.34 の三択の直後**にあり、内容は
「global CFT への最短路が欲しい読者は **THIRD PROOF** だけ読み、**Chapter I の残り全部**
(= §2-§4、Lubin-Tate 章まるごと)を飛ばしてよい」である。
★**ABC3 が既に建てた Lubin-Tate の山こそが Milne の言う「飛ばしてよい部分」**であり、
**この断りは SECOND PROOF を推してはいない**。

Milne は同じ結論に 3 本の道を用意し、見出しに依存先を書いている(実測、行番号は `.txt`):

* FIRST PROOF  `(LUBIN-TATE AND HASSE-ARF; I 2-4)`        —— 1669 行
* SECOND PROOF `(LUBIN-TATE AND COHOMOLOGY; I 2-3; III 1-3)` —— 1695 行
* THIRD PROOF  `(COHOMOLOGY AND HILBERT SYMBOLS; III 1-5)`   —— 1776 行

## ★★★★訂正 2: SECOND PROOF は Hasse-Arf を**回避しているが、Lean では FIRST より高い**

Hasse-Arf 回避そのものは**正しい**(Remark 4.15 が明言。分岐フィルトレーションも
上付き番号も 1.14/1.15 に現れない)。**しかし** 1.15 は入力として
**Chapter III Theorem 3.4**(コホモロジーによる局所 Artin 写像)を丸ごと呼ぶ。
見出し自身が `III 1-3` と書いている。下部構造の実測:

| 測定 | Ch.I §4(Hasse-Arf 路) | Ch.III §1-3(SECOND PROOF の入力) |
|---|---|---|
| 畳み込み参照 | **7 件**(硬いのは Herbrand と Hasse-Arf の 2 件) | **20 件**(Hilbert 90・Herbrand 商・Tate の定理・inf-res・カップ積) |
| 必要な mathlib 語彙 | 下付き分岐群・φ/ψ・上付き番号 —— **無いが述べられる** | `Br(K)` の乗法・`inv_K`・副有限のカップ積 —— ★**述べることすらできない** |

`memory/mathlib-cohomology-inventory-2026-09-05.md` の実測がそのまま効く:
`BrauerGroup K` には**群構造すら入っておらず、カップ積は 0 件**、`tateCohomology` は `[Fintype G]` 必須。

⇒ ★**1.5 頁は「Theorem 3.4 があれば以下は配管」という分岐器であって、費用は分岐先にある。**

## ★節点数(確定)

* SECOND PROOF 路 = **名前つき 12 + Theorem 3.4 の下部木 9 = 21 節点、うち 9 が着手不能**
  (「12」という見積は**名前つきの数としては当たり**。1 つが 9 節点の部分木だった)
* FIRST PROOF 路 = **10 節点、着手不能ゼロ**

⇒ **現時点の推奨は (B) FIRST PROOF**(mathlib に無いが全部「述べられる」)。

## ★★★思ったより安い発見が 3 つ

1. **Lemma 4.11 は Hasse-Arf を 1 滴も使わない**(2841 行)。Ch.I §4 の 3 補題のうち
   **高いのは Lemma 4.9 だけ**。4.10 は ANT 7.58(ABC3 が既に持つ不分岐拡大の一意性)1 件。
2. ★**「最大位数の巡回部分群は直和因子」という mathlib 不在の一般補題は要らない。**
   `⟨σ⟩ ∩ N = 1` と `|⟨σ⟩|·|N| = |G|` の**位数勘定だけ**で `IsCompl` が出る
   (アーベル群なので)。当初 `Module.Baer` から作る見積だった。
3. ★**SECOND PROOF の指数計算の右辺 `[K_{π,n}·K_m:K]` は既に在庫にある**
   (`Found/PGC/RamifiedUnramifiedDisjoint.lean:86` の `exists_finrank_sup_lubinTate_unramified`)。

## ★★★★★訂正 3: **Rosen 1981 は要再検討**(前回「使えない」と結論したのは早かった)

Remark 4.15 は Rosen 1981 (Trans. AMS 265) を
**「char 0 の局所体について、Hasse-Arf もコホモロジーも使わない、ただし上記より複雑な証明」**
と名指ししている。★**`PAdicLocalField p` は char 0 なので適用範囲がちょうど一致する。**
`ResearchPaper/0_Source/` に無いので**読めていない**。

- **判断待ち**: Rosen 1981 を `0_Source/` に入れるか(★人が入手する必要がある)。
  ★**Λ7 の最短化に最も効く一手**という段取り係の評価。ただし**未読なので確認できていない**。

## ★循環の注意

Example 4.7(2787 行)が「局所 Kronecker-Weber ⟹ Hasse-Arf」を示している。
つまり **Λ7 の Lean 証明は Hasse-Arf の内容を必ず含む**。含まずに済むのは
**機械(分岐フィルトレーション)だけ**。Rosen が価値を持つのはここ。

## ★原文に印の無い飛躍を 1 つ見つけた(N8)

Theorem 4.8 は Lemma 4.11(**局所体上**の主張)を**無限次拡大 `K_π` 上に適用**している。
有限段 `K_{π,n}` への降下が要る。★合図の語を持たないので `hedge-index` では出ない型。

## ★次の実装ノード(材料は全部 ABC3 の在庫内、合計 280 行の見積)

* `Found/PGC/AbelianFrobeniusSplit.lean` —— `isCompl_zpowers_frobLift`(80 行)
* `Found/PGC/AbelianSplitUnramified.lean` —— `exists_totallyRamified_sup_unramified`(200 行)

★**退化の自己検査の根拠が原典にある**: Milne **Example 4.13**(`ℚ₅` 上の具体例)が
`m` を指数より小さく取ると結論が崩れることを実証している。
⇒ `Check/PGC/Lemma411Degenerate.lean` は**原文が反例を書いてくれている珍しい 1 本**になる。

- **決定**: —

### D22・D23・D24 の合流(第 1055、★本体セッションの算術。測定値ではない)

3 つの報告が同じ数字に収束した:

* D24 の続報 ——「§2 の退化を**唯一確実に切る**のは `abelianImage`(原文の Herbrand: 像 = U_K^v)だけ」
* D23 の続報 —— Λ7 の FIRST PROOF 路のノード表で **N3 `ramificationGroup` 下付き 400 行 /
  N4 `herbrandPhi/Psi` + Prop 4.4 500 行**
* `abelianImage` を述べるには N3(下付き分岐群)が、その像が `U_K^v` だと言うには N4(Herbrand)が要る

⇒ ★**§2 の非退化化の価格 ≒ N3 + N4 ≒ 900 行**。前の agent の引き算(2,000-5,000 行)より安い。
★ただし**これは本体の算術**であって、N3/N4 の見積自体が段取り係の見積である。着手前に測ること。

★これで **§2 の非退化化と Λ7(FIRST PROOF 路)が N3・N4 を共有する**ことが確定した。
どちらを先に始めても他方の 900 行が前払いされる。

- **決定**: —

### D23 の続報(第 1055)—— Λ6(Dwork)の節点を確定した。臨界路は **13**(見積 7 は範囲が違った)

## ★★安くなった発見 2 つ

1. ★★**Dwork の逐次近似は解析ではなく純代数で書ける。**
   `ABC3.Found.GaloisRep.isAdicComplete_valuationSubring`(`Found/GaloisRep/AdicCompleteValued.lean:124`)
   に、今日着地した `isDiscreteValuationRing_unramifiedCompletionInt`(M1)と
   `CompleteSpace`(`UniformSpace.Completion` から自動)を与えると
   **`IsAdicComplete (𝔪) 𝒪_{K̂^ur}`** が出る。⇒ Cauchy 列・ε-δ を使わず
   `IsPrecomplete.prec` + `IsHausdorff.haus` で済む。M1 の `exists_eq_uniformizer_mul` が
   ちょうど「1 段進める」道具になる。
   ★**Milne の証明路(逆極限 A.7/A.8)はそのまま使えない**(ノルム完備化と 𝔪-進完備化の
   一致が未証明 —— Λ5 の逸脱が未閉)。代わりにこの道が開いた。
2. ★**乗法版 Dwork は加法版の 1 段補題を再利用できる**(`1+π^{n+1}B` の段で
   `σe−e ≡ −d` が加法版と同じ式になる)。Milne の "similar" の中身を確認した結果。

## ★Λ6a′ の正体 —— 「不要」ではなく「**Dwork から独立で、しかも安い**」

Λ6a′ = 同じ `π` の `f, g ∈ F_π` で `K_{π,f} = K_{π,g}`。これは **Prop 3.10 の `u = 1` の場合**で、
`ε = 1` が取れるので `σθ = θ`、すなわち `θ ∈ 𝒪_K[[T]]` = `[1]_{g,f}`。
**Dwork の補題を一切通らない。** 材料は全部在庫(`LubinTateEndo` / `powerSeries_uniqueness` /
`lubinTateEvalAtTorsionPoint`)。★**M3(ℂ_K)も要らない。**
⇒ 見積 400-800 → **200-400**。★**今すぐ並行で配れる独立の葉**。

## ★★★危険信号: Λ5 と Λ6 で Frobenius が食い違っている

`unramifiedClosureGalEquivZHat`(`UnramifiedZhat.lean:482`)は `coherentFrobenius`
(**位相的生成元一般**)を `Classical.choose` で使っている。Dwork は**算術 Frobenius**を使う。
`Ẑ` の同定が `Ẑ^×` のぶん不定なので、★**Λ7/Λ8 で Artin 写像を組むときに必ず衝突する。**
Λ9(捩れ)は `Ẑ` が捩れ自由なので影響なし。
⇒ 新ノード **Λ5b′ `ArithFrobeniusIsTopGen`(150-300 行)**で潰すのが安い。

## ★★`ArtinMapPiIndependent` の量化の向き(退化の予防)

`∃ π ∃ ϖ` で書くと `π = ϖ` で**自明化する**。`∀ π ∀ ϖ` でなければならない。
`LubinTateFieldPiIndependent` も同型。★**D13(`∀` が強すぎ)と逆向きの退化**であり、
退化検査ファイルを 1 本立てる価値がある(現在 12 本)。
同じく `DworkLemmaAdditive` は `∃ σ, ∀ c, ∃ b` であって `∀ c, ∃ σ, ∃ b` ではない。

## ★Theorem 3.9 の「体」の半分は Λ7 から無料で出る

Milne p.58(Example 4.13 の直後)が明言している ——
「`K_π·K^un` が `π` に依らないことは **Prop 3.10 を使わずに**回復できる。
ただし `φ` が `π` に依らないことを示すには依然としてこの命題が要る」。
⇒ **Dwork が迂回不能なのは「写像 `φ_π`」の半分だけ**(#10 `ArtinMapPiIndependent`)。
#9 `LubinTateFieldPiIndependent` は Λ7(Thm 4.8)に振り替えられる。

## ★確定した節点数: 臨界路上 **13**(★2026-09-06 に段取り係自身が訂正)

★**「見積 7 は外れ(過小)」は言い過ぎだった** —— 段取り係が背景の grep 完了後に自己訂正した。
`:1195` の 7 の内訳は `UnramifiedResidueAlgClosed`(Λ5b、済) + `DworkLemmaAdditive` +
`DworkLemmaMultiplicative` + `DworkTheta` + `DworkThetaConjugates` +
`LubinTateFieldPiIndependent` + `ArtinMapPiIndependent` であり、
★**同じ範囲を数え直しても 7 である。数え間違いではなく境界の引き方だった。**

* M1・M2・M3 は `:1059` に「Λ6 に足りない材料 3 つ」として**節点の外**に置かれていた
* Λ6a′ は `:1077` の「いま配れる持ち場」の表にあり、これも外
* ★**本当に見えていなかったのは `SubfieldClosed`(Lemma 3.12) と
  `ArithFrobeniusIsTopGen`(Λ5b′) の 2 つだけ**

CLAUDE.md の「合図 1 つ = 節点 1 つ」に従えば材料も持ち場も節点なので、
**臨界路上の総数 13** は変わらない。合計 **4,000-7,000 行**。
#9 を Λ7 に振り替えると Λ6 単独の臨界路は 11。

漏れの理由: M2 `DworkFixedRing`・M3 `ClosureCompletion`・Lemma 3.12 `SubfieldClosed`・
Λ5b′ が数えられていなかった。#9 を Λ7 に振り替えると **Λ6 単独の臨界路は 11**。
合計 **4,000-7,000 行**(D23 の 3,000-5,000 をやや上回る)。

☆段取り係が**自信が無いと明記した**箇所: M3 の 600-1,200 は幅が大きく
「`ℂ_K` 全体」と「`K̂^ur(λ)` だけ」のどちらが安いか**測っていない**。
Lemma 3.12 の 300-600 は**根拠の弱い数字**。
「位相的生成元でも `σ−1` が全射か」は**未確認**(#12 を入れれば回避できる)。
`IsAdicComplete` instance の側条件が自動で付くかも**未確認**(#1 の最初の関門)。

## ★次に配れる実装ノード 4 本(いずれも材料は在庫内)

| 優先 | ノード | 見積 | 理由 |
|---|---|---|---|
| 1 | **Λ5b′ `ArithFrobeniusIsTopGen`** | 150-300 | ★**着地済みのファイル 2 本の不整合**を潰す。放置すると Λ7/Λ8 で衝突 |
| 2 | **Λ6a′ `LubinTateFieldFIndependent`** | 200-400 | Dwork から独立の葉。M3 不要 |
| 3 | **N1+N2(Milne Lemma 4.11)** | 280 | Λ7 の FIRST PROOF 路の起点。Hasse-Arf 非依存 |
| 4 | **`DworkLemmaAdditive`** | 250-400 | Λ6 の起点。`IsAdicComplete` 経由 |

- **決定**: —

### D23 の続報 2(第 1055)—— ★`ArtinMapPiIndependent` の**必要性は「未確認」に落とす**

段取り係の自己訂正 2 点目。`:1175` に既に記録がある:

> **Λ8 には Milne は原典にならない** —— 本人が p.34 に地の文で
> 「Need to add a proof of this to the notes」と書き、鍵の可換図式を
> **Iwasawa 1986, Thm 6.9 p.89** に外注している。

段取り係は「Λ8 が要求するから #10 `ArtinMapPiIndependent`(500-1,000 行)は迂回不能」と書いたが、
★**その Λ8 側の原典が未確定である以上、#10 の必要性は「未確認」に落とすのが正直**である。

Milne p.58 の「`K_π·K^un` の π 非依存は Prop 3.10 抜きで回復できるが、
`φ` の π 非依存には依然 Prop 3.10 が要る」は **Milne 内で確かめた事実**だが、
それは「**Milne の `φ` を作るなら**」という条件つきの主張である。

⇒ ★**Λ6 の最も高いノード(#10、500-1,000 行)が本当に要るかは、Λ8 の原典が決まるまで分からない。**
Iwasawa 1986 は `ResearchPaper/0_Source/` に無い。

- **判断待ち**: Iwasawa 1986(Local Class Field Theory)を `0_Source/` に入れるか。
  ★Rosen 1981 と合わせて**入手の判断が 2 件**溜まった。どちらも人が取ってくる必要がある。

★**最初の 1 ノードの推奨は変わらない** —— `DworkLemmaAdditive` は
上の不確実性(M3 の幅・#10 の必要性・Λ8 の原典)に**一切依存しない**。
並行で配れる独立の葉は `LubinTateFieldFIndependent`(Λ6a′)。

- **決定**: —

### D22・D23・D24 の合流の**撤回**(第 1055)—— ★「§2 ≒ 900 行」は撤回する

★**本体セッションが書いた「§2 の非退化化 ≒ N3 + N4 ≒ 900 行」は撤回する。**
根拠にした N3(400 行)・N4(500 行)は、Λ7 の段取り係が背景の grep 完了後に
**「根拠のない数字だった」と自己撤回した**ものである。

★**優先すべきは D23 の実測**(`:1029`):

> **Λ7 = 10,000-18,000 行 / 8-15 ファイル**(旧見積 1,200-2,000 の 1 桁上)。
> 根拠: `RamificationGroup.lean` に `TODO: Define higher ramification groups` があり、
> 下付き `G_i` / Herbrand `φψ` / 上付き `G^v` / Herbrand の定理 / Hasse-Arf /
> 不変写像が**全部不在**。

D22 の実測では **Λ7a(下付き `G_i`)単体で 800-1,500 行**。
⇒ **§2 の非退化化は「900 行」より高い。** 正しい額は Λ7a を含む範囲で、まだ測り切れていない。

★**教訓**: 本体が別々の agent の数字を足し算するとき、**その数字の根拠の強さが揃っているか**を
見ていなかった。D23 の 10,000-18,000 は「mathlib に何が無いか」の実測に基づき、
撤回された 400/500 は agent の見立てだった。★**根拠の強さが違う数字を足さない。**

## ★Λ7 の段取り係が D23 に**足した**もの 3 点(こちらは有効)

1. ★**D23 の「Λ7 は在庫がゼロ」は言い過ぎ**。Milne **Lemma 4.11 + 4.10** は
   分岐フィルトレーションを一切使わず、材料が全部 ABC3 在庫内で **約 280 行**
   (`RamifiedUnramifiedDisjoint.lean` 110 行 + `UnramifiedSubextension.lean` 232 行の実測が基準)。
   ★**これを先に落とすと Λ7 の未知は Lemma 4.9 一本に縮む。**
2. ★**D23 の未測定欄が 1 つ埋まった** —— `:1086` の「コホモロジー経路が Milne I.4 より
   安いかは測れていない」に対し、**測った結果は「高い」**。⇒ **Milne I.4 路の選択は正しい。**
3. ★**Rosen 1981 の適用範囲が確定**。`:1083` の「実体が確認できなかった、Λ7 の最大の不確実性」に対し、
   Remark 4.15 は **「char 0 の局所体について」**と限定しており、`PAdicLocalField p` は char 0。
   ⇒ **この限定は我々には無コスト**。入手の価値は D23 の見立てより**上がった**。

## ★配れる持ち場が 2 つに割れた(独立なので同時に配れる)

| 候補 | 行数 | 消費先 | 新しい土台 |
|---|---|---|---|
| **N1+N2**(Milne Lemma 4.11+4.10) | **280**(★根拠あり) | Λ7 のみ | **不要** |
| **Λ7a**(下付き `G_i`) | 800-1,500(D23 実測) | Λ7 の Lemma 4.9 **+ §2 の 3 定理** | 要る |

★N1+N2 を先に置く理由: **Λ7 の残量が Lemma 4.9 だけになってはじめて、
`LocalReciprocityLaw` を `Interface` に置く判断が「壁の移動」でなくなる。**

## ★下流は完全に配線済み(Λ4∘Λ5 は 1 本になっている)

`Found/PGC/LubinTateZhat.lean`(56 行)の `exists_lubinTateUnramified_decomposition_zhat` が
`Gal(K_π·K^ur/K) ≃* 𝒪_K^× × Ẑ` を sorry 無しで出しており、
→ Λ8 → `ProfiniteUnitsTorsion.lean` の Λ9 まで繋がっている。
★**残っているのは Λ7 と Λ8 だけ**である。

## ★Λ6 と Λ7 の依存は「路に依存する」

D23 の `:1038`「**Λ6 は Λ7 の前提ではない**」は、選んだ路(Milne I §4)では**正しい**
—— Lemma 4.9/4.10/4.11 は `π` 非依存性に一言も触れない。
ただし SECOND PROOF の 1.14 は Theorem 3.9 を明示的に使う。
★**Ch.III 路を捨てた以上、D23 の記述はそのまま有効。**

- **決定**: —

### D24 の続報(第 1055)—— ★★§2 の価格を測った。**900 行は両方向に外れていた**

★問いが 2 つ混ざっていたのが原因。分けると額がまったく違う。

## (α) `abelianImage` を `Interface` に足して `⊤`/`⊥` を切る = **400-1,000 行、新しい数学ゼロ**

★**`Gv` は構造体が与えるデータであって構成物ではないので、N3 も N4 も 1 行も要らない。**
⇒ 本体が足した「N3 + N4 = 900 行」は**丸ごと不要**だった。

| 項目 | 見積 | 根拠 |
|---|---|---|
| 連続全射 `Γ_K ↠ 𝒪_K^×` の無条件梱包 | 80-200 | `AbsGalUnitsSurjective.lean` が 59 行で同型を組んでいる |
| `principalUnits` の橋 4 本(π 非依存・単調・`≠⊤`・`≠⊥`) | 60-150 | 既存補題の言い換え |
| `Interface` の場(`recip`/`recip_cont`/`recip_surj`/`abelianImage`) | 15-25 | statement のみ |
| witness `ltPreimageFiltration` | 120-300 | `RamificationNaturality.lean` 全体が 105 行 |
| 退化検査(13 例目) | 100-200 | `Theorem42NaiveGC.lean` が 208 行 |
| 巻き添え(`RamificationFiltration` を**構成する箇所は 3 つだけ**) | 30-100 | `topFiltration`/`topFiltration'`/`degenerateRF`。参照は 11 ファイル |

## (β) `RamificationFiltration` を本物に差し替える = **3,400-8,500 行**(原文どおり `Γ_K^{ab}` なら +Λ7)

⇒ (β) では 900 は **4-9 倍の過小**。

## ★★★★★危険信号: `abelianImage` を足すと `theorem_4_2` は「偽」から「空虚に真」へ移るだけ

`abelianImage` の witness `ltPreimageFiltration K v := comap (recip K) (principalUnits K π ⌈v⌉₊)` は
**π に本質的に依存する**(`ker recip_π = Gal(K̄/K_π)` を含み、`K_π` は π ごとに違う)。
体の自己同型は素元を動かすので、★**この witness が自然になることは原理的にない**
(Λ6 = Dwork を入れても直らない —— Λ6 が救うのは写像 `φ_π` の π 非依存性であって核ではない)。
逆に `IsNaturalFiltration` 側の既存 witness `topFiltration` は
`abelianImage` を足した瞬間に死ぬ(`map recip ⊤ = ⊤ ≠ U^2`)。

⇒ `theorem_4_2` は `(RF, hnat)` の**対**を要求するので、
★**仮説の連言に共通 witness が無い**という**13 例目の新種**になりうる。
**偽 → 空虚は前進だが解決ではない。**
★D24 の「留保 (1)」は正しかった —— 測ってみると留保どころか**ほぼ確実に起きる**。

## ★★`prop_2_2` は `abelianImage` で 1 ミリも直らない

`RF` の消費のされ方を実測した:

| 定理 | `RF` の使われ方 | `abelianImage` を足すと |
|---|---|---|
| `prop_2_2`(`Section2.lean:134`) | **`_RF` は結論にも型にも現れない** | ★**何も変わらない** |
| `cor_3_1` / `cor_3_3` | 仮説 `α : FilteredGroup.Iso` 経由 | 効くが、自由関数(`isHodgeTate`/`toGal`)の偽は残る |
| `theorem_4_2` | 結論の型 | 効く(上の空虚化の危険つき) |

⇒ ★**費用対効果が最も高いのは D24 第 1 段(自由な項関数の `∀` を外す)**。
`abelianImage` では 1 本も直らない `prop_2_2` を含む 3 本に効く。

## ★★N4 の見積は 4-10 倍の過小だった(500 → 2,200-5,400)

主因は **Prop 4.4 が Milne では `PROOF. See Serre 1962, IV.3, Pptn 14.` の 1 行**で済んでいること。
内訳: `φ`/`ψ` の区分線型構成 500-1,200 / Prop 4.4 本体 1,000-2,500 /
無限次への逆極限 300-800 / LT 塔の Herbrand(`U^{(i)}/U^{(n)} ≅ G_{q^i−1}`)400-900。
★**Hasse-Arf はこの経路には現れない**(要るのは Λ7 = Thm 4.8 の側)。
★**ここは ABC3 に地の利がある** —— `[u]_f(λ_n) − λ_n` の付値計算は既存の LT 捩れ点機構で書ける。

較正根拠(実測): `Found/PGC/` = **33,296 行 / 123 ファイル**、うち `LubinTate*` = 13,260 行 / 55 ファイル。
これが Milne CFT Ch I §2-3(約 12 頁)に相当 ⇒ **約 1,000 行/頁**。

## ★★★6 件目の「不在」の誤りを見つけて直した

`prop_2_2.needs` は「下付き番号付けも 0 件」と記録していたが、**名前が違うだけで mathlib に在る**:
`Ideal.inertia` / `AddSubgroup.inertia` / `mem_inertia`(@[simp]) /
★`subgroupOf_inertia` は **rfl**(= 原文の「下付きは部分群への移行と両立する」がそのまま在る)。
⇒ `G_i(L/K) := Ideal.inertia (L ≃ₐ[K] L) (𝔪_L^(i+1))` で**定義は 1 行**。
N3 狭義は **150-400 行**(400 という見積は妥当、やや過大)。
★`.needs` の該当箇所を**追記で訂正した**(消していない)。memory に
`mathlib-absent-by-wrong-name` を新設。

## ★`.needs` に写っていない依存を 1 つ見つけた

原文 Prop 2.1 は **`U_K ⊆ Γ_K^{ab}`(Artin 埋め込み)**を地の前提にしている
("if we regard `U_K` … as a subgroup of `Γ_K^{ab}`")。`.needs` には Herbrand と
番号付け変換しか挙がっていない。★**在庫は全射 `Γ_K ↠ 𝒪_K^×` のみで、埋め込みは Λ7 の先。**
⇒ 原文どおりの `abelianImage` は Λ7 無しには述べられない。述べられるのは
`map recip_K (Gv K v) = U^{⌈v⌉}` という**真に弱い版**(原文の帰結なので逸脱として記録すれば忠実側)。

## ☆道具の不具合を 1 件(未修正)

`hedge-index --cite` が Prop 2.1 の合図を「引用なし」と出すが、原文には
`(Theorem 1 of [3], p. 155)` がある。`.txt` で**合図が 168 行目、引用が 169 行目**にあり、
**行単位の引用抽出が行を跨げていない**。★`--cite` は「手順書」として使う道具なので影響がある。

## ★注意: `v ≤ 0` の扱い

`U^{⌈v⌉₊}` は `v ≤ 0` で `⊤` になり `Gv K v = ⊤` を強制する。
`Gv K 0 = inertia` と衝突して構造体が空になるので、★**`0 < v` に限って課すこと**。

☆段取り係が**分からなかったと明記した**こと: `⋂_π comap recip_π(U^n)` が本物の `Γ_K^v` に
一致するか(一致すれば自然な witness になる)。`Gal(K_π·K^ur/K) ≅ 𝒪_K^× × Ẑ` を使って
不分岐成分を 1 に固定できるはずだが**測っていない**。Λ7 を経由せず `Γ_K^{ab}` の
literal な形に届く抜け道は**見つからなかった**(「存在しない」の証明ではない)。

- **決定**: —

### D23 の続報 3(第 1056)—— ★★**Rosen 路は採らない。Milne I §4 を続ける**

★**本体が「Rosen が当たれば Serre は不要になるかもしれない」と述べたのは楽観的すぎた。訂正する。**

Rosen 7 頁を逐語で読んだ結果、**3 つの理由で Milne 路が勝つ**:

## ★★理由 1: 「Herbrand」が**2 つの別物**を指していた(最大の罠)

| 語 | 中身 | 出所 |
|---|---|---|
| Rosen の Herbrand | **巡回群の Tate コホモロジーの商** `h = \|Ĥ⁰\|/\|Ĥ¹\|` | Lang ANT p.188 Lemma 4 |
| [pGC] §2 / D22・D24 の Herbrand | **分岐理論の `φ`/`ψ`** | Milne Prop 4.4 = Serre 1962 IV.3 Pptn 14 |

★**両者は名前が同じだけで無関係。** Rosen は 7 頁を通して
**上付き番号付け・`φ`・`ψ`・跳び・Hasse-Arf を 1 語も使わない**(要旨の宣言どおり)。
⇒ ★**`abelianImage : map recip (G^v) = U_K^{(v)}` を述べる材料が Rosen の出力に含まれない。**
⇒ **§2 の Herbrand は Rosen 路でも別途要る。判定は YES。**

## ★理由 2: 桁は同じだが Milne が安い。しかも第 2 の消費者がある

| | 節点 | 行数 | 第 2 の消費者 |
|---|---|---|---|
| Milne FIRST PROOF | **10** | **8,000-15,000** | ★**あり** —— N3(150-400)+ N4(2,200-5,400)が §2 と共有 |
| Rosen 1981 | **16** | **10,900-24,500** | ★**無し** —— ρ13/ρ14/ρ6 の消費先は Λ7 だけ |

合算目標(Λ7 + §2 の非退化化 (β))で見ると差が開く:
Milne **8,000-15,000** 対 Rosen **13,250-30,300**。
★Milne Remark 4.15 の "more complicated" は、Lean では
**「節点 1.6 倍・行 1.9 倍」**と翻訳された。★**都合のよい方には出なかった。**

## ★理由 3: Hasse-Arf の「内容」は消えず、**未証明の外部**へ移っただけ

Rosen が回避しているのは**機械(分岐フィルトレーション)だけ**。内容は
**Borevich 1965 Thm 3,4**(`L/K` 巡回 p 次での `U^{(1)}_L` の `Z_p[C_p]` 加群構造)へ移る。
★**Rosen 自身がこれを証明していない**(行 387-397)。原論文は Proc. Steklov 80(1965、露語 16 頁)で
`0_Source` に無い。Rosen 自身が

> The proof uses the fact that `U^{(1)} = NU` in the unramified case and
> `|U^{(1)}/NU| = p` in the ramified case.

と書いており、これは **巡回 p 次でのノルム指数** —— Hasse-Arf の「跳びの整数性」と同じ鋭さ
(Milne Example 4.7 の同値性の言明と整合)。

## ★Rosen 側で見つかった安い発見 3 つ(在庫の訂正)

1. ★**正規基底定理は mathlib に在った** ——
   `exists_linearIndependent_algEquiv_apply_of_infinite`(`FieldTheory/Galois/NormalBasis.lean:62`)。
   ★**「不在」の 7 件目を回避した。**
2. ★**Noether 加群の全射自己準同型 ⇒ 単射も在った** ——
   `IsNoetherian.injective_of_surjective_endomorphism`(`RingTheory/Noetherian/Orzech.lean:60`)。
3. `Gal(Ω/K)` は既に在庫(`LubinTateZhat.lean:45`、sorry 0)。

## ★退化検査の 13 本目の候補(原文が反例を書いている 2 例目)

Rosen Lemma 11b(ii) の `Z_p × Z_p[𝔊]/N` が、**`[E:K] = p` のとき `U^{(1)}` が
`Z_p[G]` 自由でないことの実証**である。⇒ `Check/PGC/Krasner1DegenerateP.lean`。
(1 例目は Milne Example 4.13。)

## ★次の測定点として最も価値が高いところ(段取り係の指摘)

`|U^{(1)}/U^{(1)p^{s+1}}(Ψ)|` を **Borevich 無しで**出す抜け道。
完全列と Herbrand 商だけで位数が決まる可能性は排除できていない。
★もし決まるなら Rosen 路の額は 3,000-8,000 行下がり、**判定が覆りうる**。
☆ただし段取り係は「探したが見つからなかった(存在しないの証明ではない)」と明記している。

☆その他、段取り係が**分からなかったと明記した**こと: Lemma 10 の生成元 `e` の式が
OCR で落ちている(頁 5 を目視すること)/ `Lemma 11` が 2 つある(原論文の誤りか OCR の融合か不明)/
Borevich の実体を測れていない(ρ13 = 3,000-8,000 は**構成の見立て**で上下 2.7 倍の幅)。

- **決定**: —

### ★★★D25. [測定] 外部参照で畳まれた節点を全 50 本で数えた —— **水路は小さいが、止めている場所は前線そのもの**

- **状態**: **保留**(2026-09-06、第 1056)。判断は「どの文献を取得するか」

## 数え方(再実行できる形)

3 つの量を**混ぜずに**数えた。外部参照トークンは `[N]` / `[Key]` / `Author YEAR` で、
★**当該論文の書誌に実在する鍵だけ**を採る(これで Mochizuki の角括弧傍注 `[Thus]` `[Here]` が落ちる)。

* **S(外注 stub)** = `Proof.` から終端までが **30 語以下**で、書誌照合済みの外部参照を含む証明
* **N(証明なし外注)** = 見出しから次の見出しまでに `Proof` マーカーが**一度も無く**、外部参照がある項目
* **C(証明内引用)** = 証明ブロックの中の外部参照。`(項目, 文献)` の対で数える

## ★★都合の悪い側の結論: 水路自体は**小さい**

全 50 本 **15,259 証明ブロックのうち外注 stub は 56 件(0.37%)**。
★**FrdI・GenEll・pGC・CorrHyp・LocProP・NCBelyi・Falt1 で S は 1 件も無い。**
外注が集中しているのは **Milne の講義ノート 2 本で 27/56(48%)**、次いで MT の 18 件(N のみ)。

⇒ ★**「外部参照で畳む」水路が効くのは (a) Milne の講義ノート(Λ6/Λ7/Λ8 トラック)、
(b) MT / IUTchII のような引き写しの節、(c) statement 側の外注(pGC 3・CorrHyp 2)の 3 箇所だけ。**

★**ただし (c) は現に前線を止めている。件数の小ささと重要度は別である。**

## ★水路は論文ごとに違う(合計: 引用 6,639 / 語 7,415 / 傍注 11,515)

主水路は **引用 24 本・傍注 19 本・語 7 本**。同じ木の中で
**FrdI/GenEll/NCBelyi は傍注式、pGC/CorrHyp/LocProP は引用式、Milne は語式**。
★**1 本で較正した道具を全部に当てるのは今回も危険だった。**

## ★★★取得すべき文献(優先順位。下流は frontier.mjs の欄)

| # | 文献 | 塞いでいる節点 | 下流 |
|---|---|---|---|
| 1 | ★**Serre, "Local Class Field Theory"**(Cassels-Fröhlich 1967)= pGC `[3]` | ★**[pGC] Prop 1.1 は原文に証明が無い**(`Γ_K^ab ≅ (K^×)^∧` を `see, e.g., [3]` に外注)・**Prop 2.1**(`Theorem 1 of [3], p.155` = Herbrand)・Cor 1.3 | Section1 **21**(前線最大) |
| 2 | ★**Serre, "Corps Locaux" (GTM 67)** | pGC `[6]`×2、LocProP `[Ser2]`×4、★**MilneCFT `SERRE 1962`×7 = Λ7 の 2 大ノード** | Λ7 全体 |
| 3 | **Margulis** = CorrHyp `[Marg]`×7 | ★**[CorrHyp] Theorem 2.5**(証明なし・34 語)。★`frontier.mjs` が今まさにこれで `Section6.lean` を止めている | Section2 **14** |
| 4 | **Takeuchi (1983)** = CorrHyp `[Take]`×5 | [CorrHyp] Theorem 2.6(証明なし) | 同上 14 |
| 5 | **Faltings-Chai** = `[FC]` | [GenEll] Prop 3.4・Lemma 3.5 の証明内。`.needs` に 17 ファイル | GenEll 3,363 `.src` の中核 |
| 6 | **Silverman**(高さの理論)= GenEll `[Silv1] [Silv2]` | [GenEll] Prop 3.4 | `.needs` に 13 / 12 ファイル |
| 7 | **Iwasawa 1986** = MilneCFT ×3 | Λ8(Artin 写像の同変性) | 迂回不能と判定済み |
| 8 | **Lang, ANT (1970)** = FrdI `[Lang2]` | [FrdI] Theorem 6.4 の Chebotarev | — |
| 9 | **Lang 1966 + Baker 1975** | [FrdI] Lemma 6.5(六指数定理。mathlib 0 件) | — |
| 10 | ★**Mochizuki, "Foundations of p-adic Teichmüller Theory"** | MT `[Mzk2]`×**109**、FrdI ×3 | MT トラック全体 |

★**11-14**: Tamagawa 1997(LocProP ×8)/ Faltings, Crystalline(LocProP ×5)/
Elkies 1991(GenEll ×4)/ Hyodo・NSW・Fontaine・Belyi 原論文。

## ★★「取得ではなく**登記**」で済むものが 6 件ある(経路上・`0_Source` に実体あり)

`A Theory of Ordinary p-adic Curves`(FrdI `[Mzk1]`×3)/
`Categorical Representation of Locally Noetherian Log Schemes`(FrdI `[Mzk8]`×4)/
`Categories of Log Schemes with Archimedean Structures`(FrdI `[Mzk9]`×6)/
`Categories of Hyperbolic Riemann Surfaces`(FrdI `[Mzk13]`×2)/
`The Profinite Grothendieck Conjecture …`(LocProP `[Mzk1]`×5)/
`The Geometry of the Compactification of the Hurwitz Scheme`(LocProP `[Mzk3]`×2)。
★**`0_Source` の本文 104 本のうち papers.json 未登記が 54 本ある。**

★**Mochizuki 自身の論文にも本文の欠品が 3 本**(`(comments)` しか無い):
`Foundations of p-adic Teichmüller Theory` / `The Generalized Ordinary Moduli …`(RIMS 1051) /
`Combinatorialization of p-adic Teichmüller Theory`(RIMS 1076)。

## ★★測定器の欠陥 3 件(新規)

1. ★**見出し判定が ALLCAPS 書式を見ていない。** MilneCFT で見出し **2/327**、MilneANT **1/266**、
   MilneAV **1/255**。★**`--item` が引けない論文が 50 本中 15 本**ある
   (D23 が「引けなかったので手で数えた」と書いたのはこの穴)。
   ALLCAPS と Stacks のタグ前置を足すと Milne 3 本は 255/266/327 まで回復する。
2. ★★**`--cite` は `state === '済'` の項目を黙って飛ばす。** 実測の飛ばし率:
   **pGC 100%(4/4)・CorrHyp 100%(9/9)・GenEll 94%・FrdI 78%**。
   ★**我々の木が最も進んでいる論文でこそ `--cite` は何も出さない。**
3. **仏語・OCR 起因の測定不能が 8 本**(EGA1/EGA2/Del/Szp/GS/Falt1/Tate/Rosen81)。
   `Démonstration`/`Preuve` が抽出テキストに 0 件。★**語・傍注・引用のどの水路でも測れていない。**

☆測定係が**自信が弱いと明記した**こと: **N は上界**(見出し検出が漏れると区間が数千語に伸びる)。
**C は下界**(終端判定が保守的で 80 行超の証明は切れる)。
文献の所蔵判定は「著者姓がファイル名に無い」という粗い基準。

- **決定**: —

### ★★★★★D25 の続報(第 1058)—— **判定: Serre *Corps Locaux* は要る。原文が明言している**

## ★★逃げ道は原理的に無い —— Milne Remark 4.15 の逐語

> once **Proposition 4.4** and certain properties of the abelian extensions `K_{π,n}/K` are
> taken for granted, then the **Local Kronecker-Weber Theorem for K and the Hasse-Arf Theorem
> are essentially equivalent**.

⇒ ★**「Prop 4.4 + Example 4.7 を認めた上で Lemma 4.9 を証明する」= 「Hasse-Arf を証明する」。**
したがって **「LT 塔の特殊形だけで済ませる」経路は存在しない** —— 存在したらそれが Hasse-Arf の証明になる。
★**Λ7 の費用の下界が Hasse-Arf である**、と確定した。

Herbrand の側も同じ。Lemma 4.9 の本質は「`L/K` の分岐フィルトレーションと商 `K_π/K` のそれの比較」で、
Herbrand はまさにその比較定理。`φ` の推移律で代用しても、Serre 自身が
「essentiellement équivalente au théorème de Herbrand」と書いている。

## ★★本体の記述を 1 件訂正 —— Serre 1961 は「証明つき」ではなかった

★**「§3.1 が φ/ψ の初等的性質を証明つきで与えている」は過大だった。**
§3.1 の冒頭に `Nous allons rappeler` とあり**全体が復習**で、証明があるのは
COR 1・COR 2 と **`ψ` の性質 d′(整数点を整数へ、角括弧内 4 行)だけ**。
`φ` の a・b・c(連続・区分線型・狭義増加・凹・導関数)と `ψ` の a′b′c′ は**列挙のみ**。
⇒ N5/N6 は自前で証明することになる。`papers.json` の note を訂正した。

★**Serre 1961 は Hasse-Arf を証明してはいる**(§3.5 COR 1 to PROP 8)。
★**ただし §3.2 冒頭で「ここから §2 の仮定に戻る = 剰余体は代数閉」と宣言されており**、
証明も `π_0(V_n)`・`π_1(U_K^s)` のプロ代数群経由。**有限剰余体には使えない。**

## ★Artin p.99 では足りない

Serre 1961 が Artin [3] を挙げているのは **PROP 3(Herbrand)だけ**。
Hasse-Arf の代替として Artin を挙げている文献は無い(Milne も Serre 1962 V.7 のみ)。
⇒ **Artin を取っても N9(Hasse-Arf)は開いたまま。**

## ★★★ただし「取得すれば Λ7 が閉じる」ではない

Remark 4.15 のとおり Hasse-Arf は Λ7 のゴールと**本質的に同値**なので、取得は
「証明文が無くて着手できない」状態を「**2,500-6,000 行の形式化課題**」に変えるだけである。
★**取得は必要条件であって、安さの条件ではない。**

## ノード表 —— 14 ノード / **7,300-15,500 行**(前回の独立測定 10 ノード / 8,000-15,000 と同じ帯)

N1 ガロア作用の一般化 150-250 / N2 下付き `G_i` 250-400(★定義は mathlib の `Ideal.inertia`)/
N3 `𝒪_L = 𝒪_K[Π]` と `ord_L(σΠ−Π)` 判定 400-800 / N4 商の上界 400-700 /
**N5 `φ` の構成 500-900** / **N6 `ψ` と d′ 300-600** / N7 上付き 150-300 /
★**N8 Herbrand 1,200-2,500** / ★★**N9 Hasse-Arf 2,500-6,000**(★**最も根拠が弱い**。
Serre CL V.7 を読んでいない外挿)/ N10 無限次の上付き 300-600 /
N11 LT 塔の跳び 600-1,200(★**原典に完全な証明あり**)/ N12 Lemma 4.9 本体 250-500 /
★**N13 100-250** / ★**N14 200-500**。

## ★★原文が黙っている段が 2 つ(N13/N14)—— **既着地の実装に穴がある**

Milne の PROOF OF 4.8(`:2853-2862`):

1. `Lemma 4.11 applied to L·K_π` と書いてあるが、★**`L·K_π` は指数が無限**
   (`Gal(K_π/K) ≅ 𝒪^×`)で、4.11 は有限指数 `m` を要求する。
   有限水準 `L·K_{π,n}` に落として極限を取る書き換えが要る。
   ☆**原文からは Milne の意図を確定できなかった。**
2. `(4.9) implies that L_t ⊆ K_π` と書いてあるが、4.9 を使うには **`K_π ⊆ L_t`** が要る。
   4.11 の証明で `σ` を **`Gal(·/K_π)` の中から選ぶ**ことで得られるが、原文はこの選択を指定していない。
   ★★**第 1057 で着地した `exists_totallyRamified_abelian_split` はこの選択を保証していない**
   (`exists_unramified_frobenius_lift` は任意の Frobenius 持ち上げを返す)。
   ⇒ **N13 = 既着地ファイルの強化(100-250 行)。すぐ配れる。**

## ★★★危険信号: Lemma 4.9 は `K_π ⊆ L` を落とすと**偽**(13 例目の候補)

`K = ℚ_p`(p 奇)、`π = p`、`u` を非平方単数として **`L = ℚ_p(√(up))`** は
アーベル完全分岐だが `K_π = ℚ_p(μ_{p^∞})` に含まれない
(`ℤ_p^×` の平方商は位数 2 なので `K_π` 内の二次部分体は `ℚ_p(√(p*))` ただ 1 つ)。
⇒ ★**「すべてのアーベル完全分岐拡大は `K_π` に入る」という強めた形は偽。**
原文自身も 1.12(`:1633`)で近縁の現象を明示している。

★他の危険信号: `(G^n : G^{n+1}) ≤ q` を一様に書くと **`n = 0` で証明不能**
(`q = 2` で潰れない。場合分けが必須で、原文の `≤ q` はここを畳んでいる)/
★**`φ` を自由なパラメータで受け取ってはならない**(`φ = id` を入れると上付き = 下付きになり
Hasse-Arf が空になる。`G_i` から**構成**すること)/
★N2-N7 を作ったら **`Interface.PGC.RamificationFiltration` に実物として接続すること**
(接続しないと Theorem 4.2 側は `Gv ≡ ⊤` の退化 witness のまま残る)。

## ★測定器について —— 「Milne は語式」は半分しか当たっていない

§4 全体(`:2620-2870`)の合図は **8 件**、Lemma 4.9 の証明ブロック内は**わずか 1 件**。
実際の節点は 12-14。★**Milne の水路は「語」ではなく
「他項目の中の外部引用」(`See Serre 1962` の 2 つ)と「地の文の無標の断定」**である。
⇒ **MilneCFT は語式＋引用式の混合として数えるべき。**

☆段取り係が**分からなかったと明記した**こと: Serre CL V.7 の実際の長さを見ていない
(N9 の 2,500-6,000 は**外挿で根拠が弱い**。取得後に必ず測り直すこと)/
Artin が Hasse-Arf を扱うかは確認できていない/ Milne が `L·K_π` に 4.11 をどう適用する意図か/
Serre 1961 の式 (4)(5)(6)(7) は **OCR で読めず**標準形から補った/
副有限で `Subgroup.relindex` が使えるかを実測していない(N12 の見積が楽観の可能性)。

- **決定**: —

### ★D22 の未解決問いに答えが出た(第 1058)—— `↪ 𝓀` の添字範囲

D22(`:1331`)と D24 続報(`:1355`)が「`G_i/G_{i+1} ↪ 𝓀` の添字範囲が i ≥ 1 か i ≥ 0 か
記録からは読み取れない」と保留していた問いに、**3 つの原典が一致して答えた**:

* Serre 1961 `:990` PROP 1 と `:996-1004` COR 1 の証明
* Milne CFT `:2663-2666`

⇒ ★**`G_0/G_1 ↪ 𝓀^×`(i = 0 は乗法群)、`G_i/G_{i+1} ↪ 𝓀^+`(i ≥ 1 が加法群)。古典形どおり。**

★**帰結 1(D22 に効く)**: 「`i = 0` を含む形なら段差 witness も切れる」という案は**取れない**。
⇒ D24 続報の判定「**唯一確実に切るのは `abelianImage` だけ**」が**確定した**
(段差 witness `Gv K v := if v ≤ 0 then I_K else ⊥` は `↪ 𝓀` の 4 条件を通り抜ける)。

★**帰結 2(Λ7 に効く)**: Lemma 4.9 に **`n = 0` の場合分けが必須**。
原文の `(G^n : G^{n+1}) ≤ q` は `n = 0` では足りない(`q/(q−1)` が残り `q = 2` で潰れない)。
`(G^0:G^1) ∣ q−1` と `((G/H)^0:(G/H)^1) = q−1` の**割り切りの一致**で閉じる必要がある。
★**原文の `≤ q` という表示はここを畳んでいる**(合図の語を持たない断定)。

★独立の再測定が D23 と**同じ帯に落ちた**(12 ノード / 8,000-15,000 対 14 ノード / 7,300-15,500)。
D23 が「測られていない」と明記していた **Λ7a-c″ の独立見積**も出た: **2,550-4,750 行**
(本体の引き算 2,000-5,000 と整合)。

- **決定**: —

### ★D25 の続報(第 1060)—— N13 が埋まった。**Milne の意図は原文から確定できた**

## ★原文が向きを決めていた

pdftotext は `π` と `⊂/⊃` を落とすが、**Lemma 4.9 の証明 1 行目**が確定させた:

> `PROOF. Let G = Gal(L/K) and H = Gal(L/K_π), so that G/H = Gal(K_π/K).`

`H = Gal(L/K_π)` が `G = Gal(L/K)` の部分群である以上、★**仮定は `K_π ⊆ L`**(結論が `L = K_π`)。
⇒ 4.8 で 4.9 を `L_t` に当てるには **`K_π ⊆ L_t`**、すなわち **`F ≤ L_t` の向き**が要る。
逆向き(`L_t ⊆ K_π`)は 4.9 の**結論であって仮定ではない**。
☆前回「原文からは確定できなかった」と記録したが、**確定した**。

## ★★もう 1 つ穴が見えた —— N13b(無限段への移行)

4.8 の証明は `Lemma 4.11 applied to L·K_π` と書くが、
★**`Gal(L·K_π/K)` は有限指数を持たないので 4.11 は文字通りには適用できない。**
有限段 `F = K_{π,n} ≤ L·K_{π,n}` で 4.11 を回し、`n → ∞` の極限で `K_π ⊆ L_t` を出す段が要る。
⇒ **N13b**: `K_π = ⋃ K_{π,n}` に対し各 `n` の `L_t^{(n)}` の増大列から
`K_π ⊆ L_t`(無限次・完全分岐・アーベル)を作る。★**これが入って初めて Lemma 4.9 を当てられる。**
本ノード(N13)は**有限段だけ**を塞いだ。

## ★制約の代償を 1 件記録する(次の運用判断のため)

「既存の `Found/PGC/*.lean` を書き換えない」という指示のため、定理 B は
`exists_totallyRamified_abelian_split` の証明 **100 行を複製**することになった(差分は 2 か所だけ)。
見積 100-250 の上限に届いたのはこれが理由。
★**この制約はビルドの巻き添えと並行セッションとの衝突を避けるためのものだが、
「同じ証明の 2 本目」を作る代償がある。** 節目で統合するか、
強化が予見できる補題は最初から `∃` の外に仮定を出す形で書くか、判断が要る。

- **決定**: —

### ★★★★★D25 の続報(第 1061)—— **Sharifi で N8(Herbrand)の典拠が手に入った。3 分の 1 が埋まった**

ユーザーが提示した `https://www.math.ucla.edu/~sharifi/notes/algnum.pdf`
(UCLA、無償公開、259 頁)を実測した。鍵 `Sharifi` で登録済み。

## ★★★埋まったもの: **N8 = Herbrand の定理**

> **Theorem 6.5.25 (Herbrand's theorem).** `(G/N)^{φ_{L/E}(t)} = G_t N/N`

★**完全な証明つき**(Lemma 6.5.24 経由、`i_{L/K}` の言葉で書かれている)。
上付き番号版 **Proposition 9.6.6**(`(G/N)^s = G^s N/N`)も**証明つき**。

★これは **Milne CFT Prop 4.4 = Serre Corps Locaux IV.3 Pptn 14** と同一であり、
**Milne も Serre 1961 も外注していた**もの。⇒ ★**N8(1,200-2,500 行)の典拠がゼロから 1 になった。**

## ★埋まらないもの 2 つ(どちらも Sharifi 自身が "without proof" と明記)

1. **Theorem 9.6.10 (Hasse) = Hasse-Arf** —— `We state it without proof.`
   ⇒ ★**N9(2,500-6,000 行)の典拠は依然ゼロ。**
2. **Theorem 9.6.11**(`ρ_{L/K}(U_i(K)) = G^i`)—— `We next state, also without proof,`
   ⇒ ★★**[pGC] §2 が要求する `abelianImage`(= `[3]` Theorem 1, p.155)は埋まらない。**

★**重要な切り分け**: [pGC] §2 の「Herbrand」と Λ7 の「Herbrand」は**別の定理だった**:

| | 主張 | Sharifi |
|---|---|---|
| Λ7 の N8 | `(G/H)^v = G^v H/H`(商との両立) | ★**Thm 6.5.25、証明あり** |
| §2 の `abelianImage` | `ρ(U_i) = G^i`(Artin 写像による像) | ★**Thm 9.6.11、証明なし** |

⇒ **Serre / Iwasawa の必要性は §2 側と Hasse-Arf 側に残る。**

## ★その他、本文に証明つきで在るもの(我々のノード表と対応)

* **N2 下付き `G_i`** —— Def 6.5.2 / Lemma 6.5.5(正規性)/ 6.5.6 / 6.5.7(`G_i/G_{i+1} ↪ U_i/U_{i+1}`)
  / Cor 6.5.8(ガロア群は可解)
* **N3 `i_{L/K}` による判定** —— Def 6.5.14 / Remark 6.5.15(`σ ∈ G_i ⟺ i_{L/K}(σ) ≥ i+1`)
* **N5 `φ`** —— Def 6.5.19(積分形)/ Remarks 6.5.20(整数点の和の形・連続・区分線型・増加・凹・
  傾き `|G_k|/|G_0|`)/ Prop 6.5.16 / Cor 6.5.17 / Lemma 6.5.22 / 6.5.23
* **N11 相当** —— Prop 9.6.7(`ℚ_p(μ_{p^n})/ℚ_p` の上付き分岐群)+ Example 9.6.9(跳びは `0..n−1`)
  ★**円分塔での実例**なので、Lubin-Tate 塔(Milne Example 4.7)へは移し替えが要る
* 導手(Def 9.6.13, Remark 9.6.14, Example 9.6.15, Prop 9.6.16)と Lubin-Tate(§9.7)

## ★測り方(再現できる形)

WebFetch で PDF を取得 → `tool-results/` に保存されたものを PyMuPDF で抽出 →
番号付き見出しを正規表現で拾い、**直後 40 行に `Proof` があるか / `without proof` があるか**で分類。
★この粗い分類は**改頁で誤る**(Prop 9.6.6 を「証明なし」と誤判定した)ので、
**鍵になる定理は本文を逐語で読んで確かめること**。今回はそうした。

## ★取得判断の更新

| 文献 | 状態 |
|---|---|
| ~~Herbrand の証明~~ | ★**Sharifi で解決** |
| **Hasse-Arf の証明** | Serre Corps Locaux V.7。★**代替なし(Sharifi も Milne も Serre 1961 も外注)** |
| **`ρ(U_i) = G^i` の証明**(§2) | Serre "Local CFT"(Cassels-Fröhlich)[3] Thm 1 p.155。★**代替なし** |

★**講義ノートで代替できるかの判定基準が実証された** ——
「Herbrand と Hasse-Arf を**述べる**のではなく**証明している**」。
Chonoles(9 頁)は両方とも省略、Sharifi(259 頁)は Herbrand のみ証明。
★**分量が判定と相関している**(9 頁 → 0 本、51 頁 → 0 本、259 頁 → 1 本)。

- **決定**: —

### ★★★★★★★D25 の続報(第 1062)—— **Yoshida 2008 で Hasse-Arf の典拠が手に入った。経路 Λ の欠落が 3 → 1 になった**

ユーザーが MSE の質問(`Step in Yoshida's proof of Hasse-Arf theorem`)から辿った文献。
**Yoshida, "Local Class Field Theory via Lubin-Tate Theory"**、arXiv math/0606108、
Ann. Fac. Sci. Toulouse 17-2 (2008)、**21 頁**。オープンアクセス。鍵 `Yoshida08`。

## ★★★これ 1 本で埋まるもの

| 我々の欠落 | Yoshida | 状態 |
|---|---|---|
| ★★**N9 = Hasse-Arf** | **Theorem 6.11** | ★**完全な証明つき**(Sen [14] に従う) |
| N8 = Herbrand | **Proposition 6.9** | ★**証明つき**(Sharifi Thm 6.5.25 に次ぐ 2 本目) |
| ★★★**Λ7 のゴール `K^ab = K_π·K^ur`** | **Theorem 6.15**(Local Kronecker-Weber) | ★**証明つき** |
| ★★**Λ8**(Artin 写像の特徴づけと基底変換) | **Theorem A / Cor 5.16 / Thm 5.15** | ★**証明つき + 一意性つき** |

★★**N9 は、Milne も Serre 1961 も Sharifi も外注していた唯一の定理**である。
これで **典拠ゼロの節点が無くなった**(§2 の `abelianImage` を除く)。

Hasse-Arf の証明の鎖(すべて証明つき):
`Prop 6.2` → `Cor 6.3` → `Lemma 6.4/6.5` → **`Prop 6.6 (Sen)`** → `Cor 6.7` → `Lemma 6.8`
→ **`Prop 6.9 (Herbrand)`** → `Lemma 6.10` → **`Thm 6.11 (Hasse-Arf)`**。

★**Λ8 について**: Milne は p.34 に `Need to add a proof of this to the notes` と書いて
**Iwasawa 1986 Thm 6.9 p.89** に外注していた。Yoshida の要旨は
**`refining the arguments of Iwasawa [9]`** と書いており、**Theorem A (ii)** が
`Art_{K'}(x)|_{K^ab} = Art_K(N_{K'/K}(x))` を与える。⇒ **Iwasawa 1986 の取得は不要になった見込み。**

## ★埋まらないもの(1 つだけ残った)

**[pGC] §2 の `abelianImage`**(`Art(U_K^v) = Γ_K^v`、上付き分岐群と単数フィルトレーションの対応)は
**本文に見当たらない**(grep 0 件)。Yoshida は下付き `G_n` と `φ_G` だけで Hasse-Arf を通す。
⇒ **これだけが Serre "Local CFT"(Cassels-Fröhlich)[3] Thm 1 p.155 のまま。**

## ★★前提が軽く、しかも**我々の経路と一致している**

> The only prerequisites are Galois theory (including cyclotomic extensions, finite fields and
> infinite extensions) and some basic commutative algebra summarized in Appendix.

位相的なコンパクト性の議論を避けていると明記。★★**しかも Lubin-Tate 経由**なので、
我々の在庫(`Found/PGC/LubinTate*.lean` = **13,260 行 / 55 ファイル**)がそのまま土台になる。
★§3「Formal groups and Lubin-Tate groups」・§4「Lubin-Tate extensions and Artin maps」は
**我々が既に持っている部分**である。

## ★★取得判断の更新(3 → 1)

| 文献 | 前回 | 今回 |
|---|---|---|
| Herbrand の証明 | 未取得 | ★**Sharifi Thm 6.5.25 + Yoshida Prop 6.9** で解決 |
| **Hasse-Arf の証明** | ★典拠ゼロ | ★★**Yoshida Thm 6.11** で解決 |
| **Λ7 のゴール** | Milne(2 大ノードが外注) | ★★**Yoshida Thm 6.15** で自己完結 |
| **Λ8** | Iwasawa 1986(未取得) | ★★**Yoshida Thm A/5.15/5.16** で解決の見込み |
| **§2 の `abelianImage`** | 未取得 | ★**残る唯一の欠落** |

★**Rosen 棄却の判定は変わらない**(Rosen は Herbrand を使わないが `Z_p[C_p]` 格子の分類が要り、
Yoshida は 21 頁で自己完結している)。

## ★★★段取りの見直しが要る(次の測定点)

Λ7 のノード表(14 ノード / 7,300-15,500 行)は **Milne I §4 の構成**に基づいて作った。
★**Yoshida の構成は違う**(Sen の議論で Hasse-Arf を通し、`i(σ)` と `φ_G` だけで進む)。
⇒ **Yoshida の構成でノード表を作り直すと、行数が変わる可能性がある。**
特に:
* Yoshida は**上付き番号付けを最小限しか使わない**(`φ_G(n) ∈ ℤ` の形で Hasse-Arf を述べる)
* `Prop 6.6 (Sen)` が鍵で、これは `σ ∈ G_1` の位数 `p^m` に対する付値の評価
* ★**我々の在庫(Lubin-Tate 一式)がそのまま §3-§4 に対応する**

☆未確認: Yoshida の `Prop 6.2`(下付き分岐群の埋め込み)が mathlib の `Ideal.inertia` と
どう対応するか。Yoshida の `K_LT`(Lubin-Tate 拡大の合成)が我々の
`lubinTateClosure ⊔ unramifiedClosure` と一致するか。

- **決定**: —

### ★★★★★★★D25 の続報(第 1063)—— Yoshida の構成でノード表を作り直した。**取得すべき文献がゼロになる見込み**

## ★★★§2 の `abelianImage` の典拠も見つかった —— **残っていた唯一の欠落が消える**

原文本文に `Art(U^v_K) = Γ^v_K` という**定理は無い**(全文精読で確認)。
★**しかし `Prop 6.14` の証明が、その内容を有限段で完全に計算している**:

> `|G_n| = |ρ^{-1}_{f,m}(1+𝔭^i)| = q^{m−i}`  (`q^{i−1}−1 < n ≤ q^i−1`)

ここから `φ_G(q^i−1) = i` なので **`G^i = G_{q^i−1} = ρ^{-1}(1+𝔭^i) = Art(U^i_K)`**(有限段・完全分岐部分)。

⇒ ★**Serre "Local CFT"(Cassels-Fröhlich)[3] Thm 1 p.155 の取得は不要になる見込み。**
★★**これで「取得しないと着手できない」文献はゼロになった。**
ただし**別ノードとしては要る**(Y16、600-1,200 行、★根拠は弱い)——
有限段から `Γ_K^{ab}` へ上げるのに上付きの商両立 + LKW + 逆極限が要る。

## ★★Yoshida 路 = **16 ノード / 5,500-11,400 行**(Milne 路 14 / 7,300-15,500)

★**桁は同じ**(どちらも 10⁴)。中央値で Milne 11,400 / Yoshida 8,450、比 **0.74**。
★**「Yoshida なら 1 桁安くなる」は事実でない。**

**Yoshida 路で要らなくなる 5 ノード** —— ★**いずれも「原文が黙っていて我々が補った」節点**:
N6(`ψ` と `d′`)/ N10(無限次の上付き)/ **N13(着地済み 419 行)**/ N13b / N14。
⇒ **補完の負債がそのまま消える。**

**新たに要るもの**: Y2(Lemma 5.11 単項生成)/ ★**Y6(π 進展開)が最大の新規リスク**
(`C = {0}∪μ_{q−1}` の展開は Appendix I 送りで mathlib に一般形が無い)/
**Y7(Sen、Hasse-Arf の心臓)**/ Y15 の procyclic 段(`closure⟨σ⟩ ↠ Ẑ` かつ procyclic ⇒ `≅ Ẑ`)。

## ★★判定: **Yoshida 路を推す。ただし価格ではなく典拠が理由**

1. ★**証明文が全ノードに在る。** Milne 路の N9 は「文献が無くて着手できない」ままで**価格が測れない**。
2. ★消える 5 ノードが**すべて我々の補完**である。
3. ★**同じ 1 本が Λ8 の典拠でもある**(Thm 5.15 / Cor 5.16、一意性つき)。

★★**反対材料(都合のよい方に寄せない)**: **Λ8 まで含めると逆転しうる。**
Λ8 は §3-§4 の**相対** Lubin-Tate 理論を要求し、それは我々の
**54 ファイル / 13,260 行の一般化**(`[Fintype (ResidueField _)]` が **283 か所**)＋
`Θ^L_{π,π'}` ＋ Coleman 作用素である。★**この分を測っていない。**

## ★★「§3・§4 は既に持っている」は**半分誤りだった**(本体の前提の訂正)

Yoshida の §3・§4 は **`𝒪_{K̂}` 上の相対 Lubin-Tate 理論**(φ ねじれ `f∘F_f = F_f^φ∘f`)。
我々の在庫は `L = K`(`φ = id`)の**特殊化**である。`𝒪_{K̂}` の剰余体は `F̄_q` で有限でないので、
`[Fintype (ResidueField _)]` を要求する **283 か所 / 54 ファイル**はそのままでは効かない。

★**ただし Λ7 のゴール(LKW)だけなら `n = 1` で古典理論のまま走る。**
Thm 6.15 の証明は `v(σ) = n > 0` ならどの `n` でもよく、`n = 1` で `K^{ram}_x = K_π`(古典塔)。
⇒ ★**相対理論は Λ7 には要らず、Λ8 にだけ要る。**

## ★★★在庫との対応で分かった 2 つの大きな一致

* ★★**Yoshida `Prop 4.8`(`ψ : Ô^× ↠ Ô^×` が全射)= 我々の
  `DworkMultiplicative.lean::surjective_unramGalCompletionUnits_div_self`(725 行、本日着地)。完全一致。**
* ★**`Prop 4.4(iii)`(`Gal ≅ (𝒪/𝔭^m)^×`)= `galoisReciprocityEquiv`
  (`LubinTateReciprocityIsomorphism.lean:404`)。完全一致。**
* `K^{LT} = lubinTateClosure ⊔ unramifiedClosure` も一致(`n = 1` の下で)。

## ★★★既着地の 1,069 行のうち **≈360-420 行が未消費になる**(正直な数字)

| ファイル | 行数 | 判定 |
|---|---|---|
| `AbelianFrobeniusSplit.lean` | 250 | ★**全部生きる**(Y0/Y15 の `σ` の出所) |
| `AbelianSplitUnramified.lean` | 400 | 前半 ≈150 行は生きる。★**主定理 ≈205 行は Milne 4.11 専用で使わない** |
| `AbelianSplitOverSubfield.lean`(N13) | 419 | ★★**主定理 ≈190 行が生きない。** Yoshida は `Prop 5.4` から `K^{ram}_x ⊆ E_σ` を**定義から**得る |

★**取り繕わずに言うと: N13(第 1060)は Yoshida を先に読んでいれば作らずに済んだ。**
⇒ ★**N13 の追加強化(N13b)はここで止める。**
☆残り ≈650 行(Frobenius 持ち上げ・完全分岐判定・群論補題)は Y0/Y15 でそのまま消費される。
☆`AbelianSplitUnramified` の主定理は「procyclic ≅ Ẑ を避けたいときの逃げ道」として生き返りうる。
**捨てる判断はまだしなくてよい。**

## ★次の実装ノード: `Y1 = LowerRamificationGroup`(Y2 を同ファイルに畳む、550-1,050 行)

★**mathlib の `Ideal.inertia` で定義すると、原文が Prop 6.2 前半で証明している
正規性・π 非依存性が定義から無料になる**(逸脱として順序反転を記録)。
`AddSubgroup.subgroupOf_inertia` は `rfl` で、**Y7(Sen)の `H_n := G_n ∩ ⟨σ⟩` にそのまま当たる**。

★危険信号: `L/K` 完全分岐を落とすと `G_0 = ⊤` は**偽** / `Algebra.adjoin 𝒪_K {α} = 𝒪_L` を
**仮定として受け取ってはならない**(Y2 の内容が消える) / `n` と `𝔭^{n+1}` の添字ずれ /
★`G_n := ⊤`(全 n)は「大 n で `⊥`」だけを破る **退化検査の候補**。

## ☆段取り係が分からなかったと明記したこと

`Prop 4.4(i)`(`μ_{f,m} ≅ 𝒪/𝔭^m` の**加群**同型)が在庫にあるか**確定できなかった**
(`unitActionQuotientBijOn` は原始根への全単射で、全 `μ_{f,m}` の加群同型ではない)/
**α が `K(α)` の素元**という名前つき補題が見つからない(Eisenstein から従うが宣言が無い。★Y14 の起点)/
`Lemma 5.2(i)` の一般 `n` は在庫に無い(`DworkFixedRing` は `n = 1` のみ。★Λ8 で要る)/
★**Λ8 の見積を出していない**(`Fintype` 除去は 283 か所と測ったが行数は未測定)/
`Lemma 5.13` を mathlib で代替できるか型で確かめていない/
★`Prop 6.6 (iii)` の符号は**原文の見た目では決められず**、直後の整合から逆算して確定した。

- **決定**: —

## ★D26. [運用] CorrHyp の引き継ぎ時期 —— ★**決定: pGC 完了後**(2026-09-06、ユーザー判断)

## 背景(第 1063 の実測)

ユーザーから「pGC は進み始めた気がするので、ABC3b で進めていた CorrHyp を引き継げるか」と問われ、測った。

**引き継ぎ自体は技術的に可能**:

| | |
|---|---|
| CorrHyp の最終コミット | **2026-09-05 18:06**(`df863832`)—— 約 26 時間触れていない |
| 未コミットの CorrHyp 変更 | **なし**(衝突の兆候ゼロ) |
| 規模 | 21 ファイル / **11,644 行** |

★ただし**共有作業木なので、ABC3b が CorrHyp を再開しないことの確認が要る**。

## ★★塞がっている場所と進む場所が分かれている

★**`Skeleton/CorrHyp/Section2.lean` の `sorry` は 3 本だけで、3 本とも塞がっている**:

| 定理 | 塞いでいるもの |
|---|---|
| `thm_2_5` | **Margulis** p.337 Thm 27 / p.60 Lemma 3.1.1(v)。原文は**証明なし・34 語**。★未取得 |
| `thm_2_6` | **Takeuchi 1983** Thm 2.1。★未取得 |
| `prop_2_4` | 代数群の almost-simple 分解・Weil restriction・非可換 Galois コホモロジー・四元数環。 |
| | ★**mathlib の `BrauerGroup` は 7 宣言で群構造すら無い**(本日実測) |

`frontier.mjs` は `Section6.lean` を「import に現れない依存(folklore: Theorem 2.5)」で止めており、
**Section3/4/5/6 も Section2 待ち**。⇒ ★**Skeleton 側は取得なしには 1 行も進まない。**

★**一方、ABC3b が実際に進めていたのは `Found/CorrHyp/` の構成トラックで、そちらは塞がっていない。**
`df863832` の申し送り(そのまま brief になる):

> GlueData' 成分の現状(R' レベル): J・U・V・f・t・t_inv が揃った。残るのは
> **f_mono/f_open(被覆の細分待ち)と t'・t_fac・cocycle(三重交差 = pullback の同定)**。

## ★引き継いだときに進められる 3 つ(いずれも文献不要)

1. **`GlueData'` の残り**(`f_mono`/`f_open`/`t'`/`t_fac`/`cocycle`)
2. ★**NG 13 の解消** —— 現在 `check.mjs` の NG は**全件が CorrHyp の G9**(非空虚性の対照が無い)
3. ★**D17 の修理** —— `thm_6_1` の `∀ D` 量化が反証可能と確定しているが、
   「CorrHyp は不可触」を理由に**報告のみに留めていた**

- **決定**: ★**pGC 完了後に引き継ぐ**(2026-09-06、ユーザー判断)。
  理由: 実装枠は `lake` の規約で同時 1 体なので、CorrHyp を並行に走らせると pGC の波が止まる。
  ★**それまで CorrHyp は不可触のまま**(`Found/CorrHyp/**`・`Skeleton/CorrHyp/**` を触らない)。
  ★引き継ぎ時に**ABC3b の停止確認**を行うこと。

### ★★★★★★D25 の続報(第 1068)—— Y6 を測った。**原典 Lemma 6.5 は字面のままだと偽**

## ★★★★Yoshida `Lemma 6.5` の erratum(1 だけずれている)

`Def 6.1` は `i(σ) = v(σπ − π)`。原文の構成 `α_n := ∏_{i<n} σ^i(π)` から出るのは

    v(σα_n − α_n) = v(α_n) + v(σ(α_n)/α_n − 1) = n + ( i(σ^n) − 1 )

であって **`n + i(σ^n)` ではない**(原文はこの 1 行で `v(σ^n(π)/π − 1) = i(σ^n)` と読んでいるが、
正しくは `i(σ^n) − 1`)。

★**`n = 1` で決定的**: `α = π` が唯一の候補型で `v(σπ − π) = i(σ)`。しかも `v(α) = 1` なる
**任意の** `α = uπ` について `v(σα − α) = i(σ)` ちょうどになる。
⇒ ★**`σ ∈ G₁ \ {id}` では `v(σα − α) = 1 + i(σ)` を満たす `α` は存在しない。**

★★**Skeleton に逐語で書くと「偽の statement」が入る。**
この木の退化検査 12 本は「落とすと偽か自明」の検査だが、★**今回は原文の側が既に偽**という**別種**である。

★**ただし後続は壊れない。** `Prop 6.6(iii)` の消費は全体が一様に `−1` シフトするだけで、
Claim も最終評価もそのまま通る(段取り係が手で追って確認)。

★**書き方の推奨**: ℕ∞ の引き算を避け、**環の恒等式から掛け算形で出す**:
`π·(σ•α_n − α_n) = α_n·(σ^n•π − π)` ⇒ **`1 + v(σα_n − α_n) = n + i(σ^n)`**。
★`⊤` の場合分けが一切要らない(`σ = 1` でも両辺 `⊤`)。

## ★★π 進展開は「存在」だけ要る —— 一意性・収束・完備性・Teichmüller は**要らない**

| 原文が使うもの | Λ7 で要るか | 根拠 |
|---|---|---|
| 展開の**存在** | ★**要る**(回避不能) | `Prop 6.6(iii)` が項別に分解して付値の相異性を使う |
| 展開の**一意性** | ★**要らない** | `Lemma 5.2(i)` だけが使う(= Λ8) |
| **無限和・収束・完備性** | ★**要らない** | 消費側は有限の付値を否定するだけ。`M := v(z)+1` で打ち切れる |
| **`σ` の項別作用(連続性)** | ★**要らない** | 有限和なら加法性だけ(★原文は暗黙に連続性を使っている) |
| **`C = {0} ∪ μ_{q−1}`(Teichmüller)** | ★**要らない** | 完全分岐 ⇒ `𝒪_K → 𝒪_{K'}/𝔭` 全射なので代表を `𝒪_K` から取れば `σ` 不変は `smul_algebraMap` で無料 |

★**落とせないのは「零または単元」だけ**(桁が `𝔪_K` に入ると項の付値がずれ、消費側が壊れる)。10 行程度。

## ★★★負の付値も回避できる —— **分数体に ℤ 値付値を載せなくてよい**

原文は `n ∈ ℤ`(`Prop 6.6(iii)` で `s = i_{j−1} − i_j < 0` に当てる)。
`σ` が固定する元を掛けて **`s' := s + p^j·t·e` にずらすと `v_p(s') = j−1` が保たれ**、
`i((σ^p)^{s'}) = i_j` も不変。Claim も同じ論法で通る。
⇒ ★**`x, y, z, 桁`すべて `𝒪_{K'}` に収まり、`addVal : 𝒪 → ℕ∞` だけで足りる。**

★★**配管上これは大きい** —— `Found/PGC/AdjoinIntegers.lean` の冒頭が
「`IntermediateField.adjoin` に `Valued` を入れると位相のダイヤモンドで詰まった」と記録しており、
**分数体に ℤ 値付値を載せる路線は既知の地雷**である。
★付随して原文のもう 1 つの穴(`Prop 6.6(ii)` は `a ≥ 1` でしか述べていないのに `(iii)` は `a < 0` で使う)も塞がる。

★★★**この「シフト法」を Y7 の brief に必ず書くこと。落とすと鎖が切れる。**

## ★見積を根拠のある数字に置き換えた: **420-680 行**(旧 500-1,200 の下端付近)

錨: `LowerRamificationGroup.lean` = 698 行 / 上位宣言 25 本 ≒ **28 行/宣言**、
その中の `exists_mem_adjoin_sub_mem_pow` = **24 行**(Y6b の原型)。
積算 17-18 宣言で 435-515、宣言あたり 28 行での独立検算が 584。摩擦 ±30% を見て **420-680**。

⇒ ★**「Y6 が Yoshida 路の最大の新規リスク」という記述は、測った結果として支持できない。**
★**Λ7 の最大の未知は Y7(Sen (iii))に移る**(原文が 1 頁弱を割く唯一の箇所)。
☆**ただし Y7 の見積は今回測っていない。**

## ★★Y1 が既に Y6 の最大の入力を用意していた

* `exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin`
  (`LowerRamificationGroup.lean:451`、26 行)—— ★**代表系そのもの**
* `exists_mem_adjoin_sub_mem_pow`(同 `:168`、24 行)—— ★**打ち切り展開の原型**
* 作用インスタンス 3 本(同 `:536/560/582`)—— 具体層が薄くなる

★mathlib に**一般 DVR の桁展開は無い**(`PadicInt.appr` は `ℤ_[p]` 専用、
`Perfection.teichmuller` は Perfection、`WittVector.teichmuller` は Witt 環。
`representatives`/`digit`/`IsAdicComplete` で走査して該当なし)。
★`AddValuation` の「一意最小の有限和」も**乗法版しか無い**(30-45 行で自作)。

## ★`.needs` に写っていない依存が 4 本(下界であることの実例)

`σ` が付値を保つ / `v` の乗法性を `K'^×` で使う(★掛け算形にして回避)/
剰余体が伸びない ⇒ `μ_{q−1} ⊆ K` / ★**`σ` の無限和への項別作用(連続性)**(★有限打ち切りで消える)。
省略の合図は `Lemma 6.5` に **1 件**(`clearly`)だけで、この 2 補題に対応する。

## ☆段取り係が分からなかったと明記したこと

★**Y7 の行数は未測定**(「Y6 が安くなった分 Y7 が最大」と書いたが根拠は無い)/
★**シフト法を Lean で組んだ場合の摩擦は未検証**(`v_p` と ℕ/ℤ の往復で索引ずれの可能性。
★これは **Y7 の見積を押し上げる**方向)/ 退化検査 (D3) を具体例なしで書けるか未確認 /
`Lemma 6.4` は未測定 / ★**`Lemma 5.2(i)` 側の π 進展開(一意性つき・`K̂` 上)は測っていない。
「Y6 が安い」を Λ8 に外挿してはならない** / 証明本文は未構造化なので `≅` や大括弧の脱落が残る可能性。

- **決定**: —

### ★★★★D25 の続報(第 1069)—— **Λ7(Yoshida 路)のノード表を追跡ファイルに落とす**

★測定係の指摘: **この表は第 1063 の報告テキストの中にしか無く、`decisions-pending.md` にも
`memory/` にも書かれていなかった**(`grep` で 0 件)。以後の判断の土台なのでここに落とす。

典拠は **Yoshida 2008**(鍵 `Yoshida08`、構造化済み
`1_Structured/Local Class Field Theory via Lubin-Tate Theory/`、id 30 件)。

| # | 名前 | 主張 | 見積 | 状態 |
|---|---|---|---|---|
| Y2 | `MonogenicTotallyRamified` | 完全分岐 ⇒ `𝒪_{L'} = 𝒪_L[α]`(Lemma 5.11) | 250-500 | ★**済**(Y1 に同梱) |
| Y1 | `LowerRamificationGroup` | `G_n := (𝔪^{n+1}).inertia G`(Def 6.1) | 300-550 | ★**済 698 行**(Y2 込み) |
| Y3 | `RamificationQuotientEmbedding` | `G_0/G_1 ↪ 𝓀^×`、`G_n/G_{n+1} ↪ 𝓀^+`(Prop 6.2) | 400-700 | ★**済 745 行** |
| Y4 | `AbelianJumpDivisibility` | `G` 可換 ∧ `G_n ≠ G_{n+1}` ⇒ `e_0 ∣ n`(Cor 6.3) | 120-250 | 次 |
| Y5 | `SumOfConjugatesValuation` | `σ ∈ G_1 ⇒ v(Σ_{i<p} σ^i α) > v(α)`(Lemma 6.4) | 150-300 | |
| Y6 | `UniformizerExpansion` | π 進展開 + `1 + v(σα_n − α_n) = n + i(σ^n)`(Lemma 6.5) | ★**420-680**(再測) | |
| **Y7** | `SenValuationCongruence` | `i_{j−1} ≡ i_j (mod p^j)`(Prop 6.6、Sen) | 700-1,500 ★根拠弱い | ★**最大の未知** |
| Y8 | `CyclicPJumpStructure` | `G ≅ ℤ/p^m` の跳びは部分和(Cor 6.7) | 200-400 | |
| Y9 | `RamificationTransitivity` | `i(σ̄) = |H|^{-1} Σ_{τ∈H} i(στ)`(Lemma 6.8) | 400-800 | |
| Y10 | `HerbrandTheorem` | `G_n H/H = (G/H)_{φ_H(n)}`(Prop 6.9) | 250-500 | |
| Y11 | `PhiCompositionLaw` | `φ_G` の和公式と合成則(Lemma 6.10) | 400-800 | |
| **Y12** | `HasseArf` | `G` 可換 ∧ `G_n ≠ G_{n+1}` ⇒ `φ_G(n) ∈ ℤ`(Thm 6.11) | 400-900 | |
| Y13 | `UpperNumbering` | `G^m := G_{φ^{-1}(m)}`、商両立(Def 6.12 / Cor 6.13) | 350-700 | |
| Y14 | `LubinTateUpperJumps` | LT 塔の上付き分岐群(Prop 6.14) | 500-1,000 | |
| **Y15** | `LocalKroneckerWeber` | `K^ab = K_π·K^ur`(Thm 6.15)= **Λ7 のゴール** | 400-900 | |
| Y0 | `ArtinUniformizerElement` | `σ|_{K^ur}=Frob ∧ σ|_{K_π}=id`(Prop 5.4 の n=1) | 150-350 | |
| **Y16** | `ArtinUpperFiltrationImage` | `Art(U^v) = Γ^v` = ★**§2 の `abelianImage`** | 600-1,200 ★根拠弱い | |

★**着地済み 3 ノード(Y1・Y2・Y3)で 1,443 行。** 見積の合計 950-1,750 に対し実測 1,443 で**帯の中**。

### ★★Y3 で分かったこと —— **本体の brief の指定が誤りだった**

本体は Y3 の brief で「`θ_n : G_n/G_{n+1} ↪ 𝓀^+`(`(σα−α)/α^{n+1}`)が**素元の取り方に依らない**」を
求めたが、★**これは `n ≥ 1` では成り立たない。**
取り替え `β = αw` で **`c̄^α = w̄^n · c̄^β` とねじれる**(実装係が定理として証明:
`residue_ramCoeff_change_uniformizer`)。原因は **`𝔭^n/𝔭^{n+1} ≅ 𝔽_q` という同一視自身が `π^n` の選択を含む**こと。

★**原文の正規化(`𝔭^n/𝔭^{n+1}` 行き)は厳密に π 非依存**である(`thetaYoshida_independent`)。
⇒ 実装係は**両方を作り、ずれを定理化した**。消費側(Cor 6.3、Sen)は「単射」と「加法性」しか
使わないので影響なし。★**本体が Serre 型の正規化を指定したのが原因。**

★**`n = 0` / `n ≥ 1` の場合分けは原文どおり必要**で、Lean 側で**分岐点が 1 箇所に特定できた**:
`ramCoeff_mul` の補正因子 `(1 + α^n c_σ)^{n+1}` は `n ≥ 1` なら法 `𝔪` で消える
(`zero_pow (n ≠ 0)` が唯一の `hn` 使用箇所)。`n = 0` では消えない分がちょうど乗法群の構造になる。

★素元非依存性は「完全に消えた」Y1 と違い、**今回は半額**。消えたのは
「`σ ∈ G_n` という仮定が π に依らない」部分で、係数の比較 20 行は残った
(原典が `σ(π')/π' = (σπ/π)(σu/u)` と書いている計算そのもの)。

★**抽象化の切り分けが 5 回連続で効いた**: 抽象核 **0.05-0.52 秒** / 具体層 4 本まとめて 2.5 秒。

### ★Y4 の土台は揃った(次の実装)

`thetaMulQuot_injective`(`|G_0/G_1| ∣ q−1`)と `thetaAddQuot_injective` が Cor 6.3 の入力。
★実装係が**別ノードを 1 つ提案**: 「`G_1` は p 群 / `G_n/G_{n+1}` は指数 p」——
Prop 6.6(Sen)が `|⟨σ⟩| = p^m` を "by Proposition 6.2" で畳んでいる箇所の中身。
`ResidueField` の標数を `CharP` で取る必要があり、**別ノードにするのが安い**との判定。

- **決定**: —

### ★D26 の続報(第 1070)—— ABC3b / ABC3c は**停止・廃棄予定**(ユーザー通知)

2026-09-06、ユーザーより「ABC3b と ABC3c は停止しており廃棄予定」と通知があった。
`ListAgents` の実測でも **ABC3b は 2 セッションとも offline**。

## ★D26 の条件が 1 つ消える

D26 は「pGC 完了後に CorrHyp を引き継ぐ。★引き継ぎ時に**ABC3b の停止確認**を行うこと」と
決めていたが、★**廃棄されるなら停止確認は不要になる**。

★**あわせて「共有作業木での衝突」という制約も消える** ——
`autonomy-policy.md` が定める **explicit-path `git add` のみ / `-A`・`.` を使わない**、
**commit・push の前に `git fetch origin master` + `git merge-base --is-ancestor`** という規律は、
並行セッションとの衝突を避けるためのものだった。
☆ただし**規律自体は残す**(他のセッションが将来また立つ可能性があり、
explicit-path add は「意図しないものを載せない」という別の効能もある)。

## ★引き継ぎ時期の判断は変えない

★**pGC 完了後**のまま。理由は D26 のとおりで、
**実装枠が `lake` の規約で同時 1 体**だから、CorrHyp を並行に走らせると pGC の波が止まる。
⇒ 廃棄によって変わるのは「いつでも安全に引き継げる」ことだけで、**いつ引き継ぐかは変わらない**。

★**それまで `Found/CorrHyp/**` と `Skeleton/CorrHyp/**` は不可触**のまま。

- **決定**: ★**停止確認は不要になった。引き継ぎ時期は pGC 完了後のまま**(2026-09-06)。

### ★★★★★★D25 の続報(第 1072)—— **Y7(Sen)を測った。Λ7 の見積の信頼度が上がった**

## ★見積: **Y7 = 1,000-1,550 行(中心 1,250)**。★**旧 700-1,500 の下半分は支持できない**

分割を推奨: **Y7a(Prop 6.6 (i)(ii))370-610** / **Y7b((iii))610-970**。
2 通りで独立に積算し(段の積み上げ 980-1,580、宣言数 × 単価 1,000-1,500)一致した。

★**700 が出ない理由**: (iii) の 1 段落の背後に**原文が書いていない 4 本**が居る ——
「一意最小の和(★mathlib 欠落)」「打ち切り剰余」「合同の鎖」「巡回 p 群の部分群の分類」。

★**Λ7 全体 ≈ 6,140-11,730、中央 8,935**(旧中央 8,450 比 **+5.7%**)。
★**帯はほとんど動かない。最大の未知だった 1 ノードに根拠がついた。**

## ★★Y4 が前置ノードを 1 本消した

測定係が「新規ノード Y4b(`orderOf σ = p^k`、100-180 行)が要る」と測っていたものは、
★**そのまま Y4 に同梱で着地していた**(`exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one` /
`exists_nat_card_zpowers_eq_pow` = **Prop 6.6 冒頭の `(by Proposition 6.2)` そのもの**)。
⇒ ★**Y4b は削除。Y7 が自前で足すのは `m ≥ 1` の 5 行だけ。**
★もう 1 つ、上振れ条件 (c)(抽象層で `Finite G` が取れない)も **Y4 が型クラスの束を確立して消えた**。

★**Y4 の超過(見積 120-250 に対し 473 行)は「増加」ではなく「移転」**である ——
表に無かった p 群ノードが Y4 に入り、その分 Y7 から消えた。

## ★★★新しい erratum を 1 本発見 —— Prop 6.6 (iii) の帰納法の書き出し

原文 p.15 は逐語で `The assertion is empty when j = 0. **Let j = 1**, and assume the Inductive Hypothesis`。
`.txt` のバイト列でも 170dpi 画像でも `=`(同じ段落の `≥` とは字形が違う)。
★**数学的には `Let j ≥ 1` でなければ帰納法が動かない**(`j = 1` では直後の
「`σ^p` に IH を当てて `mod p^{j−1}`」が `mod p^0` になり内容が消える)。
☆**測定係は「誤植と判断したが、これは自分の推論である」と明記**(arXiv v1 と掲載版の異同は未確認)。

★さらに Claim 内の合同は **(iii) を `j−1` 以下の全添字で**要求する
⇒ **Lean では強帰納法、かつ σ について全称量化**して回すこと。

## ★★Lemma 6.5 の erratum を独立に検算した —— **(iii) は 1 か所も壊れない**

`v(σα_n − α_n) = n + i(σ^n) − 1` を `n = 0,1,2` で再確認。
使われるのは 2 箇所((A) `x` の構成、(B) 桁の像)だけで、**比較はすべて (A) 由来と (B) 由来の間**。
両辺が一様に `−1` されるので Claim も最終評価も不変。★**第 1068 の手計算を独立に再現した。**

★★**実装の推奨**: **`−1` を書かず `1 + v(…) = …` の形で通す**と、
erratum の追加コストが**実質 0 行**になる(`ℕ∞` の `1 + ·` の狭義単調性だけで往復でき、
`⊤` の場合分けも引き算も出ない)。

## ★★シフト法は通る。**しかも `e`(分岐指数)は要らない**

第 1068 の申し送りは `s' := s + p^j·t·e` だったが、必要なのは
**(a) `s' ≥ 1`** と **(b) `p^j ∣ (s' − s)`** の 2 つだけ。
`e` が要ると読めたのは「σ 不変な元を掛けてずらす」機構を想定したためで、
★**Lemma 6.5 前半は任意の `n ≥ 0` に対して直接 `α_n` を作る**ので機構を経由しない。
⇒ `i_j + s' = i_{j−1} + p^j·t` の形で **引き算を一度も書かない**。
⇒ ★**`e ≥ 1` という別ノードが 1 本消える。**

★**ℤ 値付値路線との比較(実測)**: シフト法 **+45-80 行** 対 ℤ 路線 **+250-500 行**。
差は 200-420 行で、しかも ℤ 路線には**既知の詰まり**(`AdjoinIntegers.lean` 冒頭の位相ダイヤモンド)がある。

★**シフト法の摩擦は brief の心配より小さい** —— 理由は「**`v_p` を関数として使わない**」。
原文の `v_p(…)` はすべて整除性の述語(`p^j ∣ a ∧ ¬p^{j+1} ∣ a`)に書き換わり、
★`padicValNat` / `Nat.factorization` を statement に一度も出さなくてよい
(`padicValNat p 0 = 0` というジャンク値の地雷を丸ごと回避)。

## ★★★★危険信号 —— **`ℕ∞` の引き算で書くと (iii) が空虚に真になる**(新種)

`ℕ∞` の引き算は切り詰める(`a − b = 0` when `a ≤ b`、`⊤ − ⊤ = 0`)。
★**(i) が `i_{j−1} < i_j` を与えるので、`(p^j : ℕ∞) ∣ (iSeq (j−1) − iSeq j)` と書くと
左辺が常に `0` になり `simp` で通る** —— ★**内容ゼロの定理**。
★この木の退化検査 12 本のどれとも型が違う**新種**である。

**対策**: `iN : ℕ → ℕ` の有限代表を経由し、**添字を `j+1` にして**
`(p:ℤ)^{j+1} ∣ ((iN (j+1) : ℤ) − iN j)` と書く。`⊤` は**仮定 `iSeq (j+1) ≠ ⊤` で切る**
(これが「∞ は任意の整数と合同」の正しい形式化)。

## ★(ii) は原文より浅い —— **(i) に依存しない**

`σ^a ∈ G_n ⟺ ⟨σ^a⟩ ≤ G_n`(`Subgroup.zpowers_le`)⟺ `⟨σ^{p^j}⟩ ≤ G_n` で、
★**(i) を一切使わずに** `i(σ^a) = i_j` が出る。原文の依存グラフより浅い。

## ★★本体の brief の前提が 1 つ誤っていた —— 抽出器の混同

brief は「`pdftotext` が `≠` を `=` と吐く / `⟨ ⟩` が完全に消える」を
**`0_Source/*.txt` にも当てはまる**と書いたが、実測すると `.txt`(**PyMuPDF 製**)では
**`̸=` も `⟨σ⟩` も `⊂` も残っている**(落ちるのは `Σ → P` と `≅ → ∼=`)。
★**構造化係が見たのは `pdftotext` 側**であって、抽出器が違った。
⇒ `memory/pdftotext-drops-negation` を**抽出器ごとの表に書き直した**。
☆それでも測定係は `pdftoppm -r 170` で p.14/p.15 を画像化して全文照合しており、**その姿勢は正しい**。

## ★★在庫を「型」で引き直したら 3 件ヒットした(本体の助言が効いた)

* `span_smul_uniformizer`(Y4)—— 「Cor 6.3 のために作られた」が**型は Lemma 6.5 の入力**
* `smul_mem_maximalIdeal_pow`(Y1)—— 「正規性用」だが型は「σ が付値を保つ」⇒ **印の無い依存 #8 が 10 行に減額**
* `mem_span_singleton_iff_pow_mul` / `pow_mul_mem_pow_succ_iff`(Y3)—— 「θ_n の核判定用」だが**型は桁抽出** = Y6 の核

★逆に**型で引いても無かった**もの: `addVal` は **Y3・Y4 に 1 回も現れない**。
Y1 は `addVal` 言語、Y3/Y4 は `ramCoeff` 言語で、**2 つの平行な言語**になっている。
⇒ ★**橋渡し補題 `i σ = n+1 ⟺ (σ ∈ G_n ∧ σ ∉ G_{n+1})` を 1 本立てることを推奨**(15-20 行)。

## ★mathlib の欠落を 1 件確定

**`AddValuation.map_sum_eq_of_lt`(一意最小の有限和)が無い。**
乗法版 `Valuation.map_sum_eq_of_lt` は在る。`AddValuation.` を名前欄で **全 44 件列挙**して確認。
`AddValuation` は `to_additive` ではなく**手書きミラー**なので真の欠落。自作 30-45 行。

## ★☆`.cache/decl-index.txt` が古かった(本体が直した)

**2026-09-06 14:44 生成で、Y1(21:31)/ Y3(22:09)/ Y4 を 1 件も含んでいなかった**。
⇒ agent が「不在」と誤判定する原因。★**これは「不在の誤判定」の 7 件目の回路**である
(前 6 件は名前・書式の問題だったが、今回は**索引が古い**)。
本体が作り直した(**24,039 宣言 / locator 4,610**、該当 0 → **75 件**)。
★**索引は生成時刻を見ること。**

☆測定係が分からなかったと明記したこと: "Let j = 1" が誤植か断定できない(版の異同未確認)/
Appendix I を読んでいない(★**Y6 が打ち切り形で出すかは未確認**。出ないと Y7b が上振れ)/
Y5 の行数は未測定 / ★**`lake` を動かしていないので mathlib 補題の型が合うかは未検証** /
`to_additive` 生成名は索引に載らないので乗法版から推定した / 具体層の行数の根拠が弱い /
索引が古かったため他ファイルに同等物がある可能性を排除できていない。

- **決定**: —

## ★★★D27. [運用] 実装 agent を同時 2 体にする実験 —— ★**手順をデータの前に固定する**(2026-09-07、ユーザー承認)

## 背景

`autonomy-policy.md` §4 の「`lake`/`lean.exe` を動かす agent は同時 1 体」は、
2026-09-06 に観測した**測定が壊れた 4 件**に基づく:
`lean_start` 590 秒 / `leanfile.mjs` 11-13 秒 → 1-4 分 /
★**`check.mjs` の NG が 16 ⇄ 13 で揺れる** / `setup.json` が長さ 0 で読まれる。

★**しかしメタ第 11 回が明記した**: 「`lake` 同時 1 体が律速かどうかは git からは切り分けられない。
切り分けには **2 体を同時に走らせて 1 ノードあたりの分を測る実験**が要る」。
⇒ **その実験は一度もやっていない。**

★**あわせて本体の抜けを 1 つ記録する**: メタ回が測ったのは **Λ7 の鎖だけ**で、
**pGC 全体の並列性は測っていなかった**。実際には Λ7 の外に**独立な鎖が 2 本**ある:

| 鎖 | 状態 |
|---|---|
| **Λ7**(分岐群 → Hasse-Arf → LKW) | 5/17 着地、進行中 |
| ★**Λ6 の消費側** | 数学は揃ったが**未着手 6 ノード**(`DworkThetaEval` / `SubfieldClosed` / `M3` / `LubinTateFieldPiIndependent` / `ArtinMapPiIndependent` / **Λ6a′**) |
| ★**Λ8** | 未着手。Yoshida §3-§4 の**相対** Lubin-Tate 理論が要る |

★**本体は「Λ7 が 6 本中 5 本の `sorry` の共通土台だから」直列化したが、
それは Λ7 を優先する理由であって、他を止める理由ではなかった。**

## ★測る前に固定する手順(★これを後から動かさない)

**対照(基準値)** —— すべて 2026-09-06 の実測:

* 1 ノードあたり **37.2 分**(Lean が着地した 15 session / 651 分)
* `lean_check` の往復: 抽象核 **0.05-0.6 秒** / 具体層 **0.3-2.5 秒**
* ゲート: `check.mjs --brief` **NG 13** / `selftest` **50/50** / `graph.mjs` sorry ノード **14**

**処理**: 実装 agent を **同時 2 体**にする。★**本当に独立なノードを選ぶ**
(人工的な設定を作らない) —— **Y6(Λ7)** と **Λ6a′(Λ6 の消費側)**。

**測る 3 つ**:

1. ★**ゲートの数字が揺れるか**(いちばん重要。制約の根拠がこれ)。
   両方が止まった後に `check.mjs --brief` と `graph.mjs` を回し、
   **NG が 13 のままか / ノード数が期待どおり増えたか**を見る。
   ★**揺れたら制約は正しかった**ことになり、worktree 分離へ切り替える。
2. **1 ノードあたりの所要分**(基準 37.2 分と比べる)。
3. **`lean_check` の往復時間**(基準と比べる。★役割定義が報告を求めるようになったので取れる)。

## ★★正直に記録しておくこと(交絡)

★**Y6 は単独で走り始めてから 2 体目を足す**ので、**Y6 の所要分は交絡している**。
⇒ ★**清潔なデータ点は 2 体目(Λ6a′)の側だけ**。Y6 の時間は参考値として扱う。
★ただし **(1) ゲートの安定性は重なりがあれば答えが出る**ので、そこは交絡しない。

## ★壊れたときの退避

worktree 分離(各 agent が独自の `.lake` を持つ)。代償は**冷ビルド** ——
メタ第 11 回が事故で測った限り **15 分では終わらず 5.1 GB**(M15 の 50 分と整合)。

- **決定**: ★**実験を行う**(2026-09-07、ユーザー承認「その様にしてください」)。

### ★★★★★★D27 の結果(第 1075)—— **同時 2 体でゲートは揺れなかった**

★手順は**データを見る前に固定**してある(上の D27)。以下はそれに沿った結果である。

## ★★測定 1(いちばん重要): ゲートの数字が揺れるか → **揺れなかった**

| 項目 | 事前登録した基準 | 実測 |
|---|---|---|
| `check.mjs --brief` | NG **13** | ★**NG 13 × 3 回とも同一** |
| `graph.mjs` ノード | +2 になるはず | ★**2,188(= 2,186 + 2)× 2 回とも同一** |
| `sorry` ノード | 14 | **14** |
| `selftest` | 50/50 | **50/50** |
| 文字化け | なし | なし |
| `lake build ABC3` | — | 成功(**6,986 ジョブ**) |

★2026-09-06 の事故は「**NG が 16 ⇄ 13 で揺れる**」だったので、
**3 回連続で同じ数字**が出たことが答えである。

## ★測定 3: `lean_check` の往復時間 → **基準内。異常な遅さは無かった**

★**清潔なデータ点は Λ6a′ の側**(Y6 は単独で走り始めたので交絡している):

| 項目 | 基準 | Λ6a′ の実測 |
|---|---|---|
| `lean_check` 抽象核 | 0.05-0.6 秒 | **0.09-0.76 秒** |
| `lean_check` 具体層 | 0.3-2.5 秒 | **0.67-2.12 秒** |
| `lean_start` | 10-25 秒(★**590 秒**の観測あり) | **12.0 秒** |
| 対象モジュールの `lake build` | 8-100 秒 | **12.2 秒**(下端) |
| 成果物の破損 | — | ★**兆候なし**(`.ilean` / `setup.json` の読み取り失敗も、同一コマンドで結果が変わる現象も無し) |

☆参考(交絡): Y6 は `lean_check` 0.13-1.92 秒で基準内だったが、`lean_start` が **180.8 秒**。
★ただし **Y6 の起動は 2 体目の投入前**なので、これは同時実行の影響とは言えない(冷起動の可能性)。

## ★測定 2: 1 ノードあたりの所要分 → **測れなかった(正直に)**

★**この実験では取れない。** 理由は 2 つ:
1. Y6 が交絡している(単独 → 同時)。
2. Λ6a′ 単独の所要分は取れるが、★**対照(37.2 分)は「オーケストレータの往復も含む 1 session」の値**で、
   agent の実作業時間とは定義が違う。メタ第 11 回が「**37.2 分/ノードの内訳は git からは
   切り分けられない**(M6 が未解決)」と明記しているとおり。
⇒ ★**速くなったかどうかは、この実験からは言えない。**

## ★★結論(寄せていない)

★**「`lake` 同時 1 体」の根拠だった「測定が壊れる」は、今回の 2 体では再現しなかった。**
⇒ **同時 2 体は安全に見える。** ただし:

* ★**2 体で 1 回**の観測にすぎない。3 体以上・重いノード同士では未検証。
* ★**速くなったかは測れていない**(上記)。メタ第 11 回の DAG 測定
  「3 体で飽和・上限 1.3-1.44 倍」は**依存の構造からの推定**であって、実測ではない。
* ★2026-09-06 の 4 件の観測が**何だったのか**は説明できていない。
  (`lean_start` 590 秒 / `leanfile.mjs` 11-13 秒 → 1-4 分 / NG 16⇄13 / `setup.json` 長さ 0)
  ★**別の原因(machine の負荷・他セッション・冷キャッシュ)だった可能性が残る。**

## ★運用の変更(提案。`autonomy-policy.md` §4 に当てるかは別途)

★**「独立な鎖が 2 本以上あるときは、実装 agent を 2 体まで許す」**とするのが実測に沿う。
条件:
* ★**本当に独立なノードを選ぶ**(今回は Λ7 と Λ6 の消費側)。人工的に並べない
* ★**ゲートは実装 agent が全員止まってから 1 回**(この規約は変えない。今回もそうした)
* ★**本体自身が 3 人目の `lake` 利用者にならない**(実験中は全体ビルドを回さなかった)
* ★**3 体以上は未検証**なので、当面 2 体を上限とする

- **決定**: ★**実験は完了。ゲートは揺れなかった。**★**運用変更(2 体まで)の採否は本体が節目で判断する。**

### ★★★★★★D24 の続報(第 1076)—— **3 定理とも現行形で偽であることを Lean で固定した。しかも原因は 1 つ**

`Check/PGC/FreeTermFunctionRefutation.lean`(415 行、`sorry` 0)。
3 本とも `#print axioms` は `[propext, Classical.choice, Quot.sound]`(`Skeleton` の `sorry` 定理を
参照せず**型を写した**ので `sorryAx` は入っていない)。

| 宣言 | 突いた自由な項関数 | 反例 |
|---|---|---|
| `not_prop_2_2_current_form` | `IntKbar` / `CompKbar` | `ZMod (if OneIsStandard K then 2 else 3)` + **自明な作用** ⇒ `ZMod 3 ≃+ ZMod 2` を要求して落ちる |
| `not_cor_3_1_current_form` | `isHodgeTate` | `isHodgeTate K V := OneIsStandard K` ⇒ 結論が `False ↔ True` |
| `not_cor_3_3_current_form` | `toGal` | 標準側 `toGalChoice` / 捻り側 `fun _ => 1`。★`_hρ` を**満たした上で**倒す |

## ★★★2026-09-05 の修理は 3 本とも効いていなかった

* `prop_2_2`: `SMul → DistribMulAction` の修理は効かない ——★**自明な作用はどんな型族にも乗る**
* `cor_3_1`: `V` 依存化の修理は効かない ——★**`V` を無視すればよい**
* `cor_3_3`: `_hρ` の修理は効かない ——★**満たした上で `toGal` を突ける**

## ★★★★判断材料: **3 つの穴は同一原因で、修理は 1 種類で足りる**

★**原因は「項関数に同型不変性が無い」こと**。`PAdicLocalField p` は**項**の型であって
同型類の型ではないので、`OneIsStandard K` で分岐する族が必ず作れる。
⇒ ★**修理は `ResidueCardinality.card_congr` と同じ形の `congr` 条件 1 種類で 3 本とも塞がる見込み。**
実装係が §6 にその境界を示す 2 補題を置いた
(`no_transport_badFamily` / `no_invariant_oneIsStandard` ——
「同型不変性を課せば反例族は存在しない」)。

★**これで D24 第 1 段の判断が具体化した** —— 3 本を個別に直すのではなく、**同じ条件を 1 つ足す**。

## ★既存の Check と重複していない(実測)

* `Prop22Degenerate` は**旧**形(公理ゼロの `SMul`)を非可換性で倒す。本件は**現行**の
  `DistribMulAction` 版を、**作用を一切使わずに**倒す
* `Cor33Degenerate` は**旧**形(`ρ`・`ρ'` 無関係)。本件は `_hρ` を満たした上で `toGal` を突く
* ★**`cor_3_1` は既存の反証が無く、本件が初**(`RefutationAttempts.lean` が
  「witness が作れない」で止まっていた)

★**`Skeleton/PGC/Section3.lean:98-110` の「K ≠ K′ の witness は現状の道具では構成できない」は
実測上 obsolete になった**(2026-09-05 に `twistedField` が着地して以来古かったことが、
これで確定した)。

## ★本体の見立てが外れた(良い方に)

本体は brief で「★**反例の構成なので抽象核は切り出せない可能性が高い**」と書いたが、
★**切り出せた**(`not_forall_family_addEquiv` / `not_forall_pred_iff` ほか 3 本、
**まとめて 0.90 秒**、分岐・付値・Galois の語彙が 1 つも出ない一般の添字型 `ι` の話)。
★**理由は「穴の正体が量化子の欠陥だったから」** —— 数学ではなく論理の形の問題なので一般化できた。
⇒ ★**抽象核の切り出しは 9 回連続で効いた。**

☆逸脱 3 件: `Skeleton` の `sorry` 定理を import せず型を写した / 原典の `Type*` を `Type` に固定
(★Lean では `¬ (∀ …)` の内側で宇宙変数を束縛できず不可避。`Type 0` を倒せば多相版も倒れる。
既存 2 本と同じ固定) / `cor_3_3` の反例で `E` を固定。

☆同時実行の異常: ★**無かった**(`lean_start` 1 回・13.5 秒、同じコマンドの結果のブレも無し)。
⇒ **D27 の 2 体運用の 2 例目**。

- **決定**: —

### ★★★★★D27 の重要な訂正(第 1077)—— **MCP REPL の基準環境は共有で、同時実行で壊れる**

★**2 例目の同時 2 体で、ビルド成果物とは**別の**壊れ方が出た。**

Y7a の実装係の報告(逐語):

> `lean_start(["ABC3.Found.PGC.UniformizerExpansion"])` は **10.2 秒**で「成功」と返ったが、
> その後の `lean_check` で `ABC3.Found.PGC.ramIndex` が `Unknown identifier`。
> `lean_status` を見ると imports が **`ABC3.Check.PGC.Prop12Degenerate, ABC3.Check.PGC.Cor33Degenerate`**
> (★**もう 1 体の agent のもの**)に**差し替わっていた**。
> ★**起動が 90 秒でなく 10 秒なのが合図だった。**

⇒ ★**D27 の結論を修正する。**

| 壊れ方 | 同時 2 体で | 判定 |
|---|---|---|
| **ビルド成果物**(`.ilean` / `setup.json`)とゲートの数字 | ★**壊れなかった**(NG 13 × 5 回、ノード数も一致) | 制約は不要だった |
| ★**MCP REPL の基準環境** | ★★**壊れた**(imports が他 agent のものに差し替わる) | ★**制約が要る** |

★**2026-09-06 の 4 件のうち `lean_start` 590 秒は、これで説明がつく可能性がある**
(基準環境の奪い合い)。★ただし残り 3 件は依然として未解明。

## ★対処(実装係が正しく回避した)

★**再起動しない**(役割定義どおり)。`node tools/leanfile.mjs` に切り替える
(olean を書かないので安全、**12 秒/往復**)。★実装係はこれで 3 往復で着地させた。

## ★運用の追加条件(§4 に足すべき)

★**同時 2 体のとき、MCP の `lean_start` / `lean_check` を使えるのは 1 体まで。**
2 体目は `node tools/leanfile.mjs`(ファイル単位、12 秒/往復)を使う。
★**見分け方**: `lean_start` が **90 秒でなく 10 秒で返ったら、他 agent の環境を掴んでいる**。
`lean_status` で imports を確認すること。

☆★**あるいは brief に「あなたは 2 体目なので `leanfile.mjs` を使え」と明示する**方が確実。
次の波でそうする。

- **決定**: ★**同時 2 体は継続する。ただし MCP REPL は 1 体まで**(2026-09-07)。

### ★★★★★★D25 の続報(第 1078)—— M3 を測った。**Λ6 消費側が 2,100-4,100 → 1,140-2,200(約 4 割減)**

## ★★★本体の brief の文面が**偽の statement を作るところだった**

本体は M3 の説明に「`𝒪_{K̂^ur(λ)}` が **adically complete**」と書いたが、
★★**`𝒪_{ℂ_K}` は `IsAdicComplete` では**ない**** ——
値群が稠密なので `𝔪² = 𝔪`、`⋂𝔪^n = 𝔪 ≠ 0` で **`IsHausdorff` が偽**。
⇒ 要るのは **`CompleteSpace` + `IsLinearTopology`** で、
これは `AdjoinIntegers.lean:225-284` の球イデアル構成がそのまま移る。
★同様に `IsDiscreteValuationRing 𝒪_{ℂ_K}` も**偽**。
★**`UnramifiedCompletionDVR.lean` の 447 行に対応物は無い** ——
★**これが旧見積 600-1,200 の上端が高すぎた主因。**

★**`IsAlgClosed ℂ_K`(Ax-Sen-Tate)も要らない。** 使うのは「完備な超距離ノルム体」だけ。
`Λ_n` はモニック多項式の根集合なので、`∏(θλ − ι r) = 0` から**整域性だけ**で所属が出る。

## ★判定 1: `ℂ_K` 全体は**要る**。「`K̂^ur(λ)` だけ」は**安くない、むしろ高い**

θ の係数環が `𝒪_{K̂^ur}`、評価点が `K.closure` にあり、★**両方を含む環が無ければ `θ(λ)` は型が付かない**。
`K̂^ur(λ)` を `AdjoinRoot` で作っても **`K.closure` へ戻る写像が作れない**。
⇒ ★**この選択肢は潰してよい**(段取り係の判定)。

★見積 **400-700(中心 520)**。錨: `UnramifiedCompletion.lean` の完備化構成 **128 行** +
`AdjoinIntegers.lean` の球イデアル/線形位相 **119 行** + 整数環/閉/完備 **84 行**の**逐語転写** +
mathlib の `Isometry.extensionHom`。

## ★判定 2: Λ6a′ の手は #7 に**半分だけ**効く

★**「回避」は効かない** —— Λ6a′ が閉じたのは `u = 1` で `θ ∈ 𝒪_K[[T]]` だったからで、
障害は「点の座標系」ではなく **係数環**。同じ回避は使えない。
★**「道具」は効く** —— `algHom_aeval_powerSeries_comm'`(始域≠終域版)が `𝒪_{ℂ_K}` にそのまま当たる。
★**Λ6a′ の 543 行は #7 の雛形そのもの**(抽象核 3 本が `A := 𝒪_{K̂^ur}` の条件を全部満たす)。

## ★★★Λ6 消費側の組み直し —— **約 4 割減**

| # | 旧 | 新 | 根拠 |
|---|---|---|---|
| 6 M3 | 600-1,200 ★幅大 | ★**400-700** | 代数閉性・DVR 不要。逐語転写 |
| 7 `DworkThetaEval` | 400-700 | **330-620**(+ 新規小 2 本 160-290) | Λ6a′ が雛形 |
| 8 `SubfieldClosed` | 300-600 ★根拠弱 | ★★**0(削除提案)** | ★**#9 を Λ7 に振ると消費者がゼロ**。#10 の証明は 3.12 を使わない |
| 9 | 300-600 | ★**30-80。ただし条件つき** | 下記 |
| 10 | 500-1,000 ★必要性未確認 | ★**380-800。必要性は「要る」に格上げ** | 下記 |
| **計** | **2,100-4,100** | ★**1,140-2,200(中心 ≈1,600)** | |

## ★★★#9 が無料であるための条件 —— **Y15 の量化**(★落とすと 900 行復活)

Milne p.58 が「`K_π·K^un` の π 非依存は Prop 3.10 抜きで回復できる」と保証している。
★**危険**: Yoshida `Thm 6.15` の証明は `Art_K` を使い、その `Art_K` は `Cor 4.9`(= #9+#10)で
f 非依存になっている。★**素朴に写すと #9 が Y15 の前提になり順序が逆転する。**
救い: `n = 1` では `K^{ram}_x = K_π`(古典塔)なので `Cor 4.9` を前提にせず回せる。

⇒ ★★**Y15 の brief に「結論を `∀ π, K^ab = K^ur·K_π` の形で書くこと(`∃π` にしない)」を必ず明記する。**
★**これを落とすと #9 の 300-600 行と、削除提案した #8 の 300-600 行が復活する。**

## ★#10 の必要性を「未確認」→「要る」に格上げ(★都合の悪い方向)

独立な 2 つの読みが一致した: (1) `Art_{σπ}(σa) = σ̃ Art_π(a) σ̃⁻¹` しか出ず `Gal(L/F)` 同変性に要る /
(2) ★**Λ8 の典拠 Yoshida `Cor 5.16(i)` の一意性証明が「for any uniformizer π of K」で走る**。
★ただし `Cor 5.16(i)` + `Prop 5.4` の**特徴づけ**を先に持てば #10 は 1 本に縮む。
⇒ ★**Y0(`ArtinUniformizerElement`、150-350 行)を #10 より先に置くと下がる**(額は未測定)。

## ★★構造化の判定 —— **MilneCFT はしない。Yoshida08 §2-§4 をする**

★`Found/PGC/` **139 ファイル中 `.src` を持つのは 11 本だけ**(実測)。
Found 層は元々 `.src` を要求しないので、MilneCFT 未構造化の実害は
「`brief.mjs` が使えない」だけ。296 頁の費用に見合わない。

★★**Yoshida08 は登記済みだが §2・§3・§4 が未着手**(id は §1=2 / §5=11 / §6=17 の計 30)。
★**その §3-§4(物理 p.4-9、6 頁)に Λ6 消費側の内容がそっくり入っている**
(`Def 3.3 Θ^L_{π,π'}` / `Prop 3.5` / `Lemma 4.6`(= #7) / `Prop 4.8`(**着地済み**) / `Cor 4.9`(= #9+#10))。
★しかも `index.html` 自身が「**§4 未着手 —— §6 が依存する**(Prop 6.14)」と書いており、
**Λ7 の Y14 も待っている**。
⇒ ★**優先: §4(p.6-9)> §3(p.4-6)> §2(p.2-4)。3 節あわせて 8 頁・25-35 id。**

## ★`.needs` に写っていない依存を 1 本

Milne は `f(α)=0 ⟹ g(θ(α))=0` を書くとき **`θ(α)` がどの環の元かを一言も述べない**。
合図語も無いので `hedge-index` にも出ない。★**M3 は「原典が名前を付け忘れた節点」そのもの。**
★あわせて**語表に無い合図**を 1 つ発見: 「The argument extends **without difficulty** to show … for all n」
(n=1 から一般 n への拡張)。指数に写っていない。

☆段取り係が分からなかったと明記したこと: Yoshida `Def 4.10` の 1 行が読めない(★§4 を構造化すれば解消)/
#10 の「`τ(θλ)=θλ` の計算」100-200 行は根拠が弱い / Y0 を先に置いた場合の低減額は未測定 /
★**`Completion K.closure` が `NormedField` になる instance 連鎖を実機で通していない**(★M3 の最初の関門)/
`PowerSeries.aeval` が離散性を要求しないことは仮定欄からの逆算。

- **決定**: —

### ★D27 の判定条件を訂正(第 1079)—— 「10 秒で返ったら異常」は**粗すぎた**

本体は第 1077 で「★**見分け方**: `lean_start` が **90 秒でなく 10 秒で返ったら**、
他 agent の環境を掴んでいる」と書いたが、★**これは誤り**である。

M3 の実装係(2026-09-07)の実測: `lean_start` は **10.6 秒**で返ったが
**imports は要求どおり**で、`#check` も自分の宣言が全部見えた。
⇒ ★**所要時間は合図にならない**(キャッシュが温まっていれば速く返るのは正常)。

★**正しい判定条件**: `lean_status` の **imports が自分の要求と一致しているか**を見る。
一致していれば正常、他 agent のものに差し替わっていれば異常。
`lean-idioms.md` #100 に記録済み(#98 の判定条件の訂正)。

★**運用は変えない** —— 同時 2 体のとき MCP REPL を使うのは 1 体まで、2 体目は `leanfile.mjs`。
★**ただし 1 体目も `lean_start` の後に `lean_status` で imports を確認すること。**

- **決定**: —

### ★★M3 が着地(第 1079)—— `ClosureCompletion.lean`(659 行、sorry 0)

★**測定係が「実機で通していない。M3 の最初の関門」とした instance 連鎖は 0.07 秒で越えた。**
`completableTopField_closure` を `letI` で置いた直後に `inferInstance` で出て、
★**`Valued` ダイヤモンドも起きなかった**(`Valued.v z = ‖z‖₊` は `rfl`)。
⇒ ★**見積上端 700 の理由(ダイヤモンドが荒れる)は発生しなかった。**

★**完了判定は通った**: `closureCompletionEval : PowerSeries 𝒪_{K̂^ur} →ₐ 𝒪_{ℂ_K}` が型検査を通り、
`evalAtTorsionPoint`(`Λ_n` の元での評価)まで到達。

★**抽象核が切り出せた**(9 回連続 → **10 回連続**)。
`normBallIdeal` / `isLinearTopology_of_norm_le_one` / `completeSpace_of_norm_le_one_iff` が
**`NormedField` + `IsUltrametricDist` だけ**で書けた(分岐・付値・Galois が 1 語も出ない、**0.36 秒**)。
★**同じ構成が `AdjoinIntegers.lean:225/265/278` に具体的なまま重複している**が、
着地済みなので書き換えず、**差し替え可能な形**にしてある。
⇒ ★段取り 3・4 の「128 + 119 + 84 行の逐語転写」は、抽象核を先に切ったことで **§4 が 20 行**になった。

★**`ValuationSubring` に揃えた**(理由: `unramifiedCompletionInt` が `ValuationSubring` なので
`ι_ur` の制限が同形で書ける / `IsLocalRing`・`ValuationRing`・`maximalIdeal` が無料 ——
`AdjoinIntegers` の素朴 `Subring` 版は `ValuationRing` を手証明していた)。

☆逸脱: `.src` を書いていない(MilneCFT 未構造化。★**嘘の `sectionId` は書いていない**)/
★**`IsAdicComplete` / `IsDiscreteValuationRing` / `IsAlgClosed` は作っていない**
(前 2 者は**偽**であることを理由つきで docstring に記録)/ `ℂ_K` と `K̂^ur` の関係は等長単射まで。

☆★**「他に依存しない葉」ではなくなった** —— `IsScalarTower` のために
`Algebra 𝒪_K 𝒪_{K̂^ur}` が在庫に無く、`baseIntHom`(`DworkFixedRing.lean:164`)を使うため
**`DworkFixedRing` を import した**。Λ6 消費側の入口という位置づけとは整合する。

☆新しく必要になったノード: `ℂ_K` の**代数閉性**が下流で本当に要るかの判定
(★本ファイルの想定「`∏(θλ − ι r) = 0` から整域性だけで所属を出す」が破れた場合のみ)。

---

## 第 1071 波 —— Λ7 の Sen が完結し、Yoshida §2–§4 が構造化された(2026-09-07)

### Y7b 着地 —— Prop 6.6 (iii)、これで Sen (i)(ii)(iii) が揃った

`Found/PGC/SenValuationCongruence.lean` **656 行、sorry 0**、`#print axioms` は
`[propext, Classical.choice, Quot.sound]` のみ(5 宣言すべて)。
`lake build ABC3.Found.PGC.SenValuationCongruence` 成功(2984 jobs, 7.9 秒)。

主定理 `ramIndex_pow_pow_congr`: `p^{j+1} ∣ i_{j+1} − i_j`(有限代表 `c d : ℕ` を取り ℤ で合同)。
`_or_top`(∞ 規約込み)と `PAdicLocalField` 具体層 2 本も同梱。

★**抽象核が 11 回連続で効いた**。§1+§2(`dvd_sub_of_forall_dvd_succ` /
`add_ne_add_of_not_dvd` / `eq_of_add_eq_add_of_chain` / `map_sum_ne_of_pairwise_ne` /
`smul_sum_range_sub_self`)は **ℤ の整除性と一般の `AddValuation` だけ**で書け、
**0.31 秒・一発**で通った。抽象層 §3 が 1.51 秒、具体層 §4 が 1.64 秒(差分 0.13 秒)。

★★**シフト法は Lean で本当に組めた。しかも安い。**
`s' + i_{j+1} = i_j + E` という**等式で `s'` を定義**したので、
**ℕ の引き算を一度も書いていない**。分岐指数 `e` も分数体も不要で、
`addVal : 𝒪 → ℕ∞` だけで足りた。段 7 は実質 25 行。
★原文 (ii)/(iii) の `a ≥ 1` 対 `a < 0` の齟齬も自動的に塞がる。

★★**見積より安かった** —— 610–970 行 → **656 行**(うち docstring 130 行なので Lean は約 500 行)。
段 1 は Y6 の `map_sum_eq_of_lt` をそのまま、段 2 は Y5 の
`addVal_add_le_addVal_smul_sub_self` が既にあり **新規 0 行**(見積 30–50 行)。
段 9(見積 110–180)は抽象核に落ちて 3 本 25 行 + 具体側 60 行。
★逆に段 11 の強帰納法だけは見積より厚い(`σ` の全称量化 + `σ^p` の位数を
`exists_orderOf_pow_eq` で取り直し + `j'+1 < m'` を `ramIndex_pow_pow_eq_top_iff` 経由で復元)。

★★**新しい退化species を実測で潰した** —— `ℕ∞` の切り詰め引き算では
`a < b ⟹ a − b = 0` なので `c ∣ a − b` が**空虚に真**になる
(`example (a b c : ℕ∞) (h : a < b) : c ∣ a - b` が `tsub_eq_zero_of_le` 2 行で閉じる、
前景・背景で 2 度確認)。⇒ **有限代表を取り出して ℤ で合同を書く**のが正しい形。
`lean-idioms.md` **#102**。あわせて **#101**(`rw [iff補題]` は `≠` ゴールに当たらない ——
`Did not find an occurrence of the pattern` を実際に踏んだ)。

★**原文 p.15 の "Let j = 1" 誤植の影響を受けない形に組めた** ——
シフト後の強帰納法が `j ≥ 0` 全域で回り、`j = 0` では `σ^p` への IH が
`mod p^0` で自明になるだけ。⇒ **erratum を回避する必要がなかった**。

★**MCP REPL は一切使っていない**(全往復 `node tools/leanfile.mjs`、9–11 秒/往復、往復 9 回)。
D27 の「2 体目は leanfile.mjs」規約が実際に機能した。異常なし。

### ★★Yoshida §2・§3・§4 が構造化された —— フォルダ計 68 件

| ファイル | 件数 |
|---|---|
| `section-2.html`(p.2–4) | 11 |
| `section-3.html`(p.4–6) | 14 |
| `section-4.html`(p.6–9) ★優先 | 13 |
| (既存 `section-1` 1 / `section-5` 11 / `section-6` 17) | |

構造化係が `check.mjs` の `normalize`/`squash`/`matchProjection`/`extractSections`/`PDF_MODES` を
**逐語で写した使い捨て検査**を書き、較正済み Xpdf 4.00 で照合 → **68 件中 NG 0**。
(★`section-6.html` の既存 17 件で先に一致を確認して**器具を較正**してから測っている。)
モード内訳 layout 29 / default 5 / raw 34。`mojibake.mjs` → ok。
★**`check.mjs` と `lake` は回していない**(3 人目の利用者にならないため)。

### ★★Definition 4.10 が確定した —— 測定係の推測は否だった

PDF p.9 を `pdftoppm -r 170` で目視した結果:

- ★**`K^m := K^ur L^m_f` が定義**である。`K^m = K̂^m ∩ K^sep` は**そこから従う等式**
  (f 非依存性の根拠)であって定義ではない。
- ★**「with L/K finite」は §4.1 と矛盾しない** —— `Definition 2.4` が
  「**E/K が有限なら L = Ê = E**」と明記しており、**有限不分岐拡大 `K_n` 自身が
  complete unramified extension** である。制限の理由は `L^m_f` を `K^sep` の中に置くため
  (`L = K̂` だと代数的でなくなる)。
- 論理の順序: (1) 定義 → (2) 完備化 `K̂L^m_f = K̂^m_f` が **Cor 4.9** で f 非依存 →
  (3) ★**Lemma 2.2(iii) `Ê ∩ K^sep = E`**(= Milne `Lemma 3.12` の役)で代数側へ降ろす →
  (4) `K^m` も f 非依存。
- `W(K^LT/K) ≅ W(K̂^LT/K)`(下は K であって K̂ ではない)も画像で確認。

### Λ6 消費側・Λ7 Y14 の原典対応が確定した

| この木のノード | Yoshida08 の id |
|---|---|
| #7 `DworkThetaEval` | ★**`lemma-4-6`**(p.8) |
| 着地済み `surjective_unramGalCompletionUnits_div_self` | ★**`prop-4-8`**(p.8) |
| #9 + #10 | ★**`cor-4-9`**(p.9) |
| Thm 6.15 (LKW) の §4 側 | `def-4-10` |

★**Prop 6.14 の §4 側依存先を p.17 の証明の逐行で特定した ——「3 つだけ」**:
`prop-4-4`(ii)(α が素元 ⇒ `i(σ)=v(σ(α)−α)`)・`prop-4-4`(iii)(ρ_{f,m} で翻訳)・
`lemma-4-3-ii`(`β ∈ µ^×_{f,m−i}`)。
★**`section-6.html` の依存表は「Prop 4.3(ii)」と書いていたが、原典で 4.3 は Lemma** ——
§6 側 2 箇所を訂正した。

### ★人を待たない判断として積むもの(本体が後で処理する)

1. ★★**`surjective_unramGalCompletionUnits_div_self` の `.src` が Milne LEMMA 3.11 を指している**が、
   Yoshida08 では `prop-4-8` である。証明の骨格(逆極限 → `θ^{q−1}` 全射 → Artin–Schreier)は
   Lean 側と一致しているので**間違いではない**が、Yoshida を主軸に据えた今は
   `prop-4-8` を指すほうが整合する。★**今は書き換えない**(既存 `Found/PGC/*.lean` を触ると
   走行中の agent の import が再ビルドになる)。**ゲート後に単独で行う。**
2. ★**`class="statement example"` を規約に足すか** —— 構造化係が `ex-3-6` `ex-3-8` `ex-4-2` で
   使ったが、README §3 の第 2 語の一覧に `example` が**無い**。既存コーパスに
   非 legacy 3 ファイル・4 件の先例がある。⇒ **README §3 に追記する**(実態が先行している)。
3. ★**`papers.json` の `Yoshida08` に `verifiedPages` フィールドがそもそも無い**。
   README §2-3 は記録を要求している。p.2–9 を今回 170dpi で全面目視したので
   **1–17 を入れる作業が残っている**。
4. `section-6.html` を 2 箇所変更した(依存表の訂正・末尾 `p.open` の「§4 未構造化」解消)。
   `.verbatim` は不変。**ゲートで S1–S6 を確認する。**

### ★★新たに実測した脱落 —— ハットは版面上の大きさで 2 通りに壊れる

★**本文サイズでは丸ごと落ちて `K`、上付き・下付きサイズでは帽子が `b` という文字になって残る**
(`Θ^{K̂,×}` → `Θ b K,×`、`𝒪_{K̂}` → `O b K`)。
★★**`Lemma 2.2(iii)` は帽子が落ちると `E ∩ K^sep = E` という自明な式に見える。**
他に `∐`/`⋃`/`∏` は出力なし、`⟼`→`−→`・`↦`→`→`(縦棒が消える)、`⟹`→`=⇒`、
**行末ハイフネーションが本文に残る**(p.4 `ex-amples`)、`µ` は U+00B5。
★**`⊂` は落ちない**(本体が brief に書いた説明と実測が食い違った点。brief の誤り 7 件目)。
`≠`→`=` は Prop 4.8 で再現し `data-txt="="` で明示した。

⇒ [[pdftotext-drops-negation]] のメモに**ハットの 2 通りの壊れ方**を足すこと。

### ★記法の復元に失敗／不完全と申告された箇所(取り繕わない)

- **太字 µ・太字 𝔭**: 原典は `\boldsymbol\mu` を一貫して使うが対象の区別に効いていないので
  装飾クラスを付けていない(README §5「同定に関わらない装飾は落としてよい」)。
- **`∐` と `⋃` の字形判定**: どちらも `pdftotext` 出力が空なので照合では区別できない。
  画像で Prop 4.7(i) は `∐`、それ以外は `⋃` と読んだが**低解像度では紛らわしい**。
- **`φ` の 2 種**: 算術 Frobenius `\varphi` の出力は **U+03D5**、§6 の Herbrand 関数 `\phi_G` は
  **U+03C6**。表示側も出力に合わせたが**字形そのものは検証していない**。
- ページ境界で原典の 1 項目を 2 単位に割った(`Lemma 4.3` → `lemma-4-3` + `lemma-4-3-ii`、
  `Definition 3.2` → `def-3-2` + `setup-3-1-hom`)。両ファイルに明記済み。
- 証明本体は未収録(§6 と同じ方針)。§2.1 のノルム付値の式(分数の横線が落ちる)も未収録。

### ★ゲートで確認すべきこと(本体の宿題)

Y7b が走行中に見た `check.mjs --brief` は **NG 42** だった(平常は NG 13)。
★**構造化係が `1_Structured/` を編集している最中に測った値**なので、
**確定した数字ではない**。★**ゲートで測り直すまで「NG が増えた」と言わないこと。**
期待: ノード 2,190 → 2,192(M3 + Y7b。走行中の #7・Y8 が着地すればさらに +2)、
sorry ノード 14(不変)、`1_Structured: S1–S6 すべて PASS`、G1 の新規 NG 0。

---

## 第 1072 波 —— Cor 6.7 が着地し、brief が原典の Proof を出すようになった(2026-09-07)

### Y8 着地 —— Corollary 6.7（跳び位置の p 進展開）

`Found/PGC/SenJumpExpansion.lean` **448 行(うち証明は約 180 行)、sorry 0**、
`#print axioms` は主要 4 本すべて `[propext, Classical.choice, Quot.sound]`。
`lake build ABC3.Found.PGC.SenJumpExpansion` 成功(2985 jobs, 7.8 秒)。往復 7 回・各 11–13 秒。

主定理 `exists_padicDigits_card_lowerRamificationGroup`:
`∃ nn : ℕ → ℕ, (∀ i, 1 ≤ nn i) ∧ ∀ j, 1 ≤ j → j+1 ≤ m → ∀ n,`
`Nat.card (G_n) = p^(m−j) ↔ (Σ_{i<j} nn i * p^i < n ∧ n ≤ Σ_{i≤j} nn i * p^i)`。

★★**本体の机上検算「`n_0 := i_0 − 1`」は正しかった**(訂正なし)。
Lean の中で `i_j = 1 + Σ_{i≤j} n_i p^i`(`j < m`)が出て、原文の `Σ_{i<j} < n ≤ Σ_{i≤j}` が
`i_{j−1} ≤ n < i_j` にちょうど化けた。★**`n_0 ≥ 1` は `σ ∈ G_1` から 1 度だけ使って出る**
—— つまり原典が `n_i ∈ ℤ_{≥1}` と書いたこと自体がこの読み方を強制していた、という読みが確認された。

★★**巡回群の位数から部分群を復元する補題は「不要になった」** —— mathlib を探す前に
別ルートが見つかった: `|G_n| = p^{m−(j+1)} < p^m ⟹ σ ∉ G_n ⟹ i_0 ≤ n`、`i_m = ⊤ ⟹ ¬(i_m ≤ n)`、
よって §1 の `exists_index_boundary` で `i_k ≤ n < i_{k+1}` なる `k < m` が取れ、
順方向を当てて `Nat.pow_right_injective` で `k = j`。
★**部分群の分類も `IsCyclic` の一意性も使っていない**(見積 40–80 行 → **実測 20 行**)。
⇒ ★**本体が brief で「ここだけ見積が外れやすい」と警告した箇所が、逆に最も安かった。**

★**抽象核が 12 回連続で効いた**。`exists_padicDigits_of_lt_of_dvd_sub`(純 ℕ/ℤ の望遠鏡和)+
`exists_index_boundary`(純 ℕ)が **1 往復 12 秒・修正 1 回**、
純群論の `orderOf_pow_pow_eq` / `exists_generator_of_isCyclic` が **1 往復 11 秒・一発**。

逸脱 6 件(docstring に記載)。主なもの: `hσ : σ ∈ G_1` は**仮定**(Prop 6.2 未形式化)/
`G ≅ ℤ/p^mℤ` を `(htop : zpowers σ = ⊤, hord : orderOf σ = p^m)` の対に翻訳/
`1 ≤ j ≤ m−1` → `1 ≤ j ∧ j+1 ≤ m`(ℕ の切り詰め引き算回避、端点は落としていない)。
`lean-idioms.md` に **#103**(`omega` は `p^(k+1)*c` と `c*p^(k+1)` を別原子として数える)・
**#104**(述語引数の補題は呼ぶ側が `(P := fun k => …)` を書く)。

### ★★`brief.mjs` に「1b. 原典の Proof 段落」を足した(+約 90 行)

**動機**: Y7a と Y7b の実装者が**独立に**「brief に Proof が無く、`0_Source/*.txt` を
直読して段落を取ったのが決定打だった」と報告した(Y7b は 981–1030 行を直読して
12 段の段取りをそのまま得た)。★**同じ報告が 2 回出たら道具にする**、という運用。

`0_Source/<file>.txt` から**行頭の見出し / `Proof` / `□` で切り出す**。
★**必ず行番号を添えて出す**(発見的な切り出しなので、読み手が直読できるようにする)。

**実測(Yoshida08 の 68 項目)**: ★**Proof が出た 31 / 原典に Proof なし 20 /
原典に見出しを持たない単位(我々が切り出した `setup`・`remark`)17**。
★Prop 6.6 では **986–1021 行**が出て、Y7b の実装者が手で探し当てた 981–1030 とほぼ同じだった。

**この機能が効く論文(実測)**: 行頭に見出しが立つ 11 本 ——
`BK CorrHyp EtTh GenEll LocProP NCBelyi pGC SemiAnbd Stacks Tate Yoshida08`。
`MilneCFT`(154 件)・`Sharifi`(973 件)は**走り込み式**なので退避経路で拾い、`runin` を立てて断る。
★**効かない 2 本**: `Falt1`(スキャン由来で見出し 0 件。`elementaryfibration` のように語が繋がる)・
`DelSB616`(仏語 `Lemme 1.6.-` かつ文字化け)。★**これは道具の限界であって直せない。**

### ★★[pGC] は `Proof:` ブロックを持たない —— 論証が主張の**手前**にある

★実測: pGC 33 項目のうち **Proof が出たのは 4 件だけ**。しかしこれは欠落ではなく、
★**Mochizuki が地の文で導いてから `Proposition 1.2:` と宣言する書き方**だからである
(`.txt` 95–114 行が Prop 1.2 の論証そのもので、"Thus, it follows from the above isomorphism that"
で終わって主張に入る)。
⇒ ★**`noproof` のときは「直前の地の文」を 40 行を上限に出す**ようにした。
★**これは我々の北極星（pGC Prop 1.1/1.2）に直接効く。**

★あわせて直した 3 件:
1. ★**見出しの正規表現が末尾のピリオドを弾いていた** —— `Definition 2.4.` の "." が
   否定先読み `(?![0-9.])` に当たって**行頭の見出しを見落とし、走り込み退避が
   「Lemma 4.6, ...」という本文中の引用を掴んでいた**。`(?!\.?[0-9])` に直した。
   ★**退避経路は「見つからない」より悪い誤りを作りうる**という実例。
2. `setup` / `remark` / `data-implicit` の単位は「見出しが無い」ではなく
   **「原典に見出しを持たない単位（欠落ではない）」**と言わせる。
3. `itemKeyOf`(check.mjs G1 と共有。**触っていない**)は数字番号しか拾わないので、
   `Theorem A` のような英字番号は `data-item` から拾い直す。

★`--json` にも `proof` を足した。`--node` モードは無変更。

### ★人を待たない判断として積むもの(前波からの継続 + 追加)

前波の 1–4 は未処理のまま。追加:

5. ★**`brief.mjs` の Proof 抽出は 200 行で打ち切る**。Yoshida08 では打ち切り 0 件だったが、
   長い証明を持つ論文では効く可能性がある。★**上限を上げるかは、実際に打ち切りが出てから決める**
   (今上げると出力が膨らむだけ)。
6. ★`GenEll` は Proof 16 件すべてが「終端記号なし」だった。
   ⇒ **論文ごとに証明終端の記号が違う**。`GenEll` の終端が何かは未測定。

---

## 第 1073 波 —— Λ6 消費側の入口が開いた + brief の境界判定を直した(2026-09-07)

### #7 着地 —— `DworkThetaEval`（= Yoshida `Lemma 4.6`、Milne Prop 3.10 の消費側）

`Found/PGC/DworkThetaEval.lean` **807 行(見積 450–840 の範囲内)、sorry 0**。
`#print axioms` は主要 7 本すべて `[propext, Classical.choice, Quot.sound]`。
`lake build ABC3.Found.PGC.DworkThetaEval` 成功(3431 jobs、当該モジュール 27 秒、警告 0)。
★**#7a / #7b / #7c すべて入った**(分割せずに済んだ)。

★**抽象核が 13 回連続で効いた**。`subst_iteratedLubinTate_of_semilinear_intertwine`
(`CommRing A` + `ψ : A →+* A` だけ)が **0.15 秒**、
`exists_mem_roots_of_eval_map_eq_zero`(`IsDomain` だけ)が **0.06 秒**、
`subst_intertwine_of_comp_inverse` が **0.07 秒**。具体層は 1.18–3.07 秒。

★★**代数閉性なしで根の同定ができた**(段取りの想定どおり)。
`exists_mem_roots_of_eval_map_eq_zero` は `IsAlgClosed` を一切要求せず、
`Q = ∏(X−r)` を写して `∏(θλ − ι r) = 0` にし、`Multiset.prod_eq_zero_iff` だけで所属を出した。
⇒ ★**M3 の報告で挙がっていた「`ℂ_K` の代数閉性が要るか」という新ノードは、不要と確定した。**

★★**半線形 `aeval` の橋は要らなかった**(段取りの「新規 40–70 行」は**不要**)。
`f, g` を `map (baseIntHom K)` で先に `𝒪_{K̂^ur}` 係数へ上げると絡みの両辺が**同じ係数環**になり、
既存の `aeval_subst_eq_aeval_aeval` がそのまま当たる。半線形性は `(PowerSeries.map ψ)^[n]` として
**冪級数の側だけに閉じ込めた**。★段取りより安い方向の差分。

★**M3 の `closureCompletionEval` はそのまま使えた**(`rfl` で通った)。
足りなかったのは 1 点だけ ——「`‖θ(λ)‖ ≤ ‖λ‖` ゆえに値がまた `HasEval`」で、新規約 25 行
(`norm_aeval_closureCompletionInt_le` / `norm_evalAt_lt_one` / `hasEval_evalAt`)。
☆**M3 に引っ越す価値がある**が、走行中の判断ではないので**動かしていない**。

★**原典の省略の合図を 1 つ潰した**: 51 頁 "The argument extends without difficulty … for all n" を、
はじめから任意の `n` で証明した。
逸脱 6 件(docstring に記載)。`.src` は書いていない(MilneCFT は `1_Structured/` に無い。
★**嘘の `sectionId` は書いていない**)。
`lean-idioms.md` に **#105**(statement に instance が要るなら `haveI` では間に合わない ——
`attribute [local instance]`)・**#106**(`Polynomial.coe_map` は無い。正しくは `Polynomial.polynomial_map_coe`)。

★同時実行の異常なし。`lean_start` 1 回(26 秒)→ `lean_status` で imports 一致を確認 →
以後 `lean_check` 38 回、`lean_reset` 0 回。最後の `lean_status` でも imports は自分の 3 本のまま。
⇒ ★**D27 の「MCP REPL は 1 体だけ・2 体目は leanfile.mjs」規約が 2 度目も機能した。**

### ★新しく必要になったノード —— `Lemma 4.5`（実測して確認した）

Cor 4.9 の証明は**ちょうど 3 つ**を消費する: Prop 4.8(θ の存在、着地済み)・
**Lemma 4.6**(体の一致、= #7 が今日着地)・**Lemma 4.5**(`[θ]^{(j)}[xπ_j] = [xπ′_j][θ]`)。
★**Lemma 4.5 は木に無い**(`decl-index.txt` 24,255 宣言に対して 0 件。
`DworkTheta.lean:42` の docstring に名前が出るだけ)。
⇒ ★**Cor 4.9 の後半(`ρ_{f,m} = ρ_{f′,m}`)は Lemma 4.5 待ちで、前半(体の一致)だけを #9 として配った。**

### ★★`brief.mjs` の境界判定を直した —— 「行頭に見出し語がある」だけでは境界にならない

★**実測した不具合**: Yoshida `Corollary 4.9` の証明の中に、折り返しのせいで
`Proposition 4.8. Lemma 4.6 shows bKm` という**行頭に見出し語が立つ行**ができ、
**見出しと誤認されて証明が 6 行のうち 2 行で切れていた**。
★道具自身は「終端記号に当たらず切った」と警告を出していたので**気づけた**
—— ☆**警告を出す設計にしておいたことが効いた実例**。

**直し方**: ★**各見出し鍵の「最初の出現」だけを境界として採る**
(項目は原典の中で昇順に一度だけ立つ、という性質を使う)。`boundary` / `isBoundary`。

**退行の測定**（件数は動かず、「終端記号なし」だけが減った）:

| 論文 | 項目 | Proof が出た | 原典に Proof なし | 見出しを持たない単位 | 終端記号なし |
|---|---|---|---|---|---|
| Yoshida08 | 68 | 31 | 20 | 17 | **2**(旧 5) |
| GenEll | 36 | 16 | 16 | 0 | **4** |
| pGC | 33 | 4 | 29 | 0 | **0**(旧 4) |

☆★**正直な留保**: 「終端記号なし」が減ったのは**切り出しが伸びた**ことを意味するが、
**伸びた先が正しいかは 1 件(cor-4-9、516–529 行・`□` で終わる)しか目視していない。**

### ★メタ第 12 回を起動した(隔離 worktree)

★**メタ係は第 11 回以来動いていなかった**ので起動した。持ち場は 3 つ:
(1) **M22・M24**(索引の正規表現が非 ASCII / `_root_` で壊れる。M14 の 3 例目・4 例目、未着手)
—— ★`.cache/decl-index.txt` は 24,255 宣言まで育ち、**全実装者の唯一の在庫入口**になっている。
(2) `brief.mjs` の新機能を**残り 10 論文**で測る(本体は 4 本しか測っていない)。
(3) ★**本体の brief の誤り 7 件のうち、機械で事前に捕まえられたのは何件か**の測定のみ。

★規約どおり `lake build` と `check.mjs --lean` を禁じた(M35: 隔離 worktree で `--lean` を叩くと
**cold build 50 分**が始まる)。★MCP REPL も禁じた。
★**本体が `brief.mjs` を worktree 作成後に変えた**ので、その旨と新しい基準値を伝えてある
——**採用は「差分の提案」の形で受ける**(ファイル全体の置き換えは採用しない)。

### Y9 着地 —— Lemma 6.8（商で測った i は剰余類上の平均）

`Found/PGC/QuotientRamIndexAverage.lean` **478 行、sorry 0**。
`#print axioms` は 5 本すべて `[propext, Classical.choice, Quot.sound]`。
`lake build ABC3.Found.PGC.QuotientRamIndexAverage` 成功(2986 jobs, 9.0 秒)。
★往復 **7 回**、失敗は **1 回だけ**。MCP REPL は使わず全往復 `leanfile.mjs`(11.99 秒/往復)。

主定理 `card_mul_ramIndex_eq_sum_ramIndex`:
`(Nat.card H : ℕ∞) * ramIndex ϖ σ = ∑ τ : H, ramIndex π' (σ * τ)`。
★★**除算も切り詰め引き算も 1 箇所も書いていない**(掛け算形。#102 の処方どおり)。

★★**商群 `G ⧸ H` を作らずに済んだ ——「一段まるごと消える」という段取りの見込みが当たった。**
`MulSemiringAction (G ⧸ H) C` は 1 行も書いていない。`G` を `C` に直接作用させ、
`hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c` を仮定に置くだけで `i(σ̄)` が書ける。
★★**副産物: `H` の正規性 `G ▷ H` も明示的な仮定として要らなかった** ——
正規性が担っていた内容(σ が `K″` を保つ)は `[MulSemiringAction G C]` と
`hcomp`(`algebraMap C B` の `G`-同変性)に吸収されている。

★**抽象核が 14 回連続で効いた**。`associated_map_sub_self_prod`(一般の可換整域 `R`、
`φ : R →+* R` のみ。分岐・付値・Galois の語彙が 1 つも出ない)が **1 往復 12 秒で一発**。
具体層・付値層も各 1 往復で一発。原文の 2 段(`a ∣ b` / `b ∣ a`)が
`dvd_prod_sub_map_of_coeff` と `prod_sub_map_dvd` の 2 本にそのまま写った。

★**見積 500–900 行 → 実測 478 行(思ったより安い)**。安くなった要因 3 つ:
(1) 商群を作らなかった (2) ★**「`g(X) − π″` が `f` で割れる」を仮定に置かずに導けた**
(`Multiset.prod_X_sub_C_dvd_iff_le_roots` + `Multiset.le_iff_subset` +
`eq_one_of_smul_eq_of_adjoin_eq_top` による根の相異性。段取りは「仮定でよい」という含みだった)
(3) `B^H = 𝒪_{K″}` と `𝒪_{K″} = 𝒪[π″]` を `hfix` 1 本に畳めた。
★段取りが挙げた在庫 7 本のうち**実際に使ったのは 3 本**だけ。

★**`addVal_prod`(`addVal (∏) = ∑ addVal`)は mathlib に無かったので Y9 が作った。**
`lean-idioms.md` **#105**(`Algebra.adjoin_singleton_eq_range_aeval` の `obtain` は
`.toRingHom` を被せて返すので `have hq' : (aeval t) q = y := hq` を 1 行挟む)。

### ★★Y9 が仮定に置いたもの = 新しいノード（Y10 として配った）

```
hram : ∀ c : C, addVal B (algebraMap C B c) = (Nat.card H : ℕ∞) * addVal C c
```
★**Y9 は在庫を両側(ABC3・mathlib)で測ったうえで「無かった。仮定に置いた。導いていない」と
正直に書いた。** ★mathlib の `Ideal.ramificationIdx` は
`IsDiscreteValuationRing.addVal` と接続されていない。

⇒ **Y10 `FixedRingRamificationIndex` として配った**。
★`.src` は書けない(Yoshida はこれを独立の主張として立てておらず "totally ramified" に畳んでいる)。
★★**本体は 2 つの道(共役の積=ノルム / Artin + e·f=n)を挙げたが、どちらが閉じるかを
実機で確かめていない**ことを持ち場に明記し、**測って選ぶこと・閉じなければ正直に報告すること**を
指示した。★段 1(「離散付値の延長は素元での値だけで決まる」)だけでも価値がある形で分割可にした。

### ★`brief.mjs` に字形の凡例を足した（Y9 の実装者の要請）

Y9 の報告: 「1b は**決定的に役に立った**。段取りを 1 分で確定でき、`lean-search` を 1 回も
呼ばずに済んだ。原文の `f^σ(π′) = ±b` の **`±`** が明示されていたので、
符号を単位として吸収する補題を**最初から**用意でき、符号で往復する事故が起きなかった」。

★**改善要求もそのまま出てきた**: 「`Q τ∈H(στ(π′)−π′)` の `Q` が `∏` だと**推測**できたから
読めた。凡例があれば推測が要らない」。
⇒ ★**実測して凡例を足した**（`Q`=∏ / `P`=Σ / `L`=⊕ / `S`=⋃ / `T`=⋂ / `b`+大文字=ハット /
`∼ =`=≅ / `̸=`=≠）。★**引用の中に実際に出ているものだけ**を表にして出す(noise を出さない)。
★**単独の大文字は本当にその文字のこともある**(Yoshida は体の名前に `L` を使う)ので、
曖昧なものは「★曖昧」と明示する。
あわせて「`.txt` には `===== [page N] =====` が行として挿入されるのでページ境界で文が割れる」
(Y9 が「by Lemma / 5.11.」で踏んだ)を注記に入れた。

☆Y9 の 3 つ目の指摘「『腑に落ちなければ `.txt` を直読』の誘導は今回不要だった」は
**採らなかった**——1b は発見的な切り出しであり、直読の逃げ道を消すのは危険なため。

### #9 着地 —— Corollary 4.9 の**前半**（体の f 非依存性）

`Found/PGC/LubinTateTowerFIndependent.lean` **722 行(うち docstring 約 270 行)、sorry 0**。
`#print axioms` は `[propext, Classical.choice, Quot.sound]`。
`lake build ABC3.Found.PGC.LubinTateTowerFIndependent` 成功(3432 jobs, 25 秒、警告 0)。

主定理 `lubinTateCompletionField_eq`:
`lubinTateCompletionField … f … n = lubinTateCompletionField … g … n`
(`:= IntermediateField.adjoin (unramifiedCompletion K) (closureCompletionCoe K '' ↑Λ_{f,n})`
= `ℂ_K` の中の `K̂^ur(Λ_{f,n})` = 原典の `K̂^m_f`)。★**両方の包含が入った**(`le_antisymm`)。

### ★★#8 `SubfieldClosed` の決着 —— **「0 行」は誤りだった。約 45 行である**

段取りでは「#8 → 0 行(削除提案)。完備体の有限次拡大は完備だから自動で出るはず」と
**推測**していた。★**#9 が実機で測った結果、0 行ではなく約 45 行**。
ただし段取りが恐れた 250 行でもない。

★★**決め手は 1 点**: `Submodule.closed_of_finiteDimensional` は **`NormedSpace` を要求しない**
(`[NontriviallyNormedField 𝕜] [CompleteSpace 𝕜] [Module 𝕜 E] [ContinuousSMul 𝕜 E]
[T2Space E] [IsTopologicalAddGroup E]` だけ)。
⇒ `NormedSpace K̂^ur ℂ_K`(`‖a•x‖ = ‖a‖‖x‖`)を組む必要がない。
★**これは名前でなく型で引いて初めて分かった**(`#check` 一発、0.15 秒)。**在庫を型で引く」の 9 度目の実証。**

自前で置いたのは 5 本だけ(`Algebra` と `IsScalarTower` の scoped instance 2 本 /
`ContinuousSMul` 3 行 / `closureCompletionAlgHom` + `isIntegral_closureCompletionCoe` /
抽象核 `isClosed_intermediateField_of_finiteDimensional` 8 行)。

★**抽象核が 15 回連続で効いた**。`aeval_mem_of_isClosed`(**閉部分環は冪級数の値を含む**、
0.28 秒)は分岐・付値・Lubin-Tate・Galois の語彙がゼロ。
★★**段取り係の提案「全単射 ⇒ adjoin ⊇」ではなく、「閉部分環 + 評価点 + 係数環の像 ⇒
冪級数の値も入る」という形に切り出された** —— 全単射では足りないという段取りの警告は正しく、
`trunc N θ` の極限で吸収する形が正解だった。
★実際 `exists_bijOn_iteratedLubinTateTorsionPoints`(#7 の完成形)は**使われなかった**。
使われたのは `exists_mem_torsionPoints_evalAt` と `coe_aeval_evalAt_comp_inverse` の 2 本。

★**Λ6 の結論には `constantCoeff θ' = 0` が含まれていない**ことが判明し、
新補題 `constantCoeff_of_subst_eq_X`(`θ∘θ′ = X` と `θ(0)=0` から)を用意して導いた。

`lean-idioms.md` に **#107**(★`open … in` は docstring の**前**に置く。
`lean_check` では docstring を省くので**ファイルに書いた瞬間に初めて出る** ⇒
**docstring を足したら `lake build` の前にもう一度 `lean_check` に通すこと**。#105 と同型)・
**#108**(「有限次元なら閉」は `NormedSpace` を要求しない。
`Algebra.fg_adjoin_of_finite` は根名前空間)。

☆逸脱 4 件。うち 3 番: **原典に無い「整数側の影」を併記した**が、
`Subring.topologicalClosure` を取った形で述べており、`𝒪_{K̂^ur}[Λ]` **自身**が閉であることは
扱っていない(原典に対応するのは体の側で、そちらは有限次元だから自動的に閉)。

### ★brief の誤りの 8 件目にならずに済んだ例

★**段取りが「#8 は 0 行かもしれないが確かめていない。測ること」と書いた**ので、
実装者は測ってから進み、45 行という正しい答えを出した。
☆**「推測を推測として渡す」ことが機能した実例**として記録する
(過去 7 件の誤りは、いずれも推測を断定として渡したものだった)。

### ★実装者からの `brief.mjs` の指摘 2 件 —— 1 件は既に直っており、1 件は**採らない**

#9 の実装者は**私が境界規則を直す前の版**を見ており、`cor-4-9` の切り出しが 4 行で
切れたと報告した。現在は **516–529 行**(`□` まで全文)が出る。

- 提案 1「行全体が `□` を終端に足す」→ ★**既に入っている**(`END` が行末の `□` に当たる)。
- 提案 2「**直前の行が空行でない見出し候補は無視する**」→ ★★**採らない。実測で否定された。**
  Yoshida の `.txt` では**見出しの直前は常に本文**である(85 行目 `Definition 2.1.` の前は
  "…E = S K′ = Frac(OE)."、88・104・129・142・179・196・218 行も同様)。
  ⇒ この規則を入れると**見出しが 1 つも取れなくなる**。
  本体が入れた「**各見出し鍵の最初の出現だけを境界にする**」規則の方が正しい。
- 提案 3「行番号を出しているのは非常に良い。この 1 点で往復が 1 回で済んだ」→ 維持する。

### 新しく配った持ち場: `Lemma 4.5`（Cor 4.9 の後半を解禁する）

★**Cor 4.9 の後半 `ρ_{f,m} = ρ_{f′,m}` は Lemma 4.5 だけを待っている。**
★抽象核は **1-コサイクル**である(原文が "argue by induction in both directions" と
書いているのがそれ) —— `a (j+j') = σ^j (a j') * a j` を満たす 2 つが `j = 1` で一致すれば
`∀ j : ℤ` で一致する、という純群論。
★**`j` は Weil 群の Frobenius 次数なので負の `j` が本当に要る**ことを明記した。
★在庫: **`uniformizerProd (σ : G) (π : B) (n : ℕ)`(`UniformizerExpansion.lean:164`)が
Yoshida の `π_m = ∏_{t<m} π^{ϕ^t}` と同じ形**で既にある(★型で引いて見つけた)。
ℤ へ延ばすか ℤ 版を新設するかは**本体が確かめていない**ので、測って選ぶよう指示した。

### ★Prop 6.9（Herbrand）の先読み —— **主張そのものより先に必要な節点が 2 つある**

原典(`.txt` 1054–1065 行、証明は 6 行で完全)を本体が直読した。

> Proposition 6.9 (Herbrand). Define **φ_H(n) := −1 + (1/|H|) Σ_{τ∈H} min{i(τ), n+1}** for **n ∈ ℝ≥0**.
> Also, for n ∈ ℝ≥0, define **G_n := {σ ∈ G | i(σ) ≥ n + 1}**。Then **G_n H/H = (G/H)_{φ_H(n)}**.

★★**この木の `lowerRamificationGroup` は ℕ 添字である**(`((𝔪_B)^(n+1)).inertia G`)。
Herbrand は **実数添字**を要求する。⇒ **先に立てる節点が 2 つ**:

1. ★**実数添字の `G_n`**(`{σ | i(σ) ≥ n+1}`、`n : ℝ≥0`)。整数 `n` で Definition 6.1 と
   一致することも要る(原典が "i.e. `G_n = G_i` if `i ∈ ℤ≥0` and `n ∈ (i−1, i]`" と明記)。
2. ★★**`φ_H` の定義 —— ここで初めて「除算」が避けられない。**
   `φ_H(n) = −1 + (1/|H|) Σ min{i(τ), n+1}` は**本質的に有理数値**である。
   ⇒ ★**`ℕ∞` から出て `ℝ`(または `ℚ`)へ移る決断がここで要る。**
   ★これまで(#102 以来)「除算を書かない・掛け算形で持つ」で通してきたが、
   **Herbrand 関数はその方針が原理的に通らない最初の場所**である。
   `i(τ) = ⊤`(τ = id)があるので `min{⊤, n+1} = n+1` の扱いも要る
   ⇒ `ℝ≥0∞` を経由するか、`i` を有限値に落としてから和を取るか。★**未決**。

★**Y9 の掛け算形はここで正しく効く**: 原典の証明は Lemma 6.8 を
`i(σ̄) = φ_H(m−1) + 1` の形で使うが、Y9 の
`(Nat.card H : ℕ∞) * ramIndex ϖ σ = ∑ τ, ramIndex π' (σ*τ)` に
`i(στ) = min{i(τ), m}` を代入すれば **`|H| · i(σ̄) = Σ min{i(τ), m}`** となり、
`φ_H` の定義式と**除算を挟まずに**突き合わせられる。
⇒ ☆**除算は `φ_H` の定義の中だけに閉じ込められる見込み**(★本体は実機で確かめていない)。

★証明が使うもう 1 つの材料: **`i(στ) = min{i(τ), m}`**。
これは `i(τ) ≥ min{i(στ), i(σ^{-1})}`(★**i の超距離不等式**)から出る。
★**この不等式が木にあるかは未測定。**

⇒ **次に配る持ち場は Herbrand 本体ではなく「実数添字の `G_n` と `φ_H`」**とする。
★Y10(`e = |H|`)が走行中で実装枠が埋まっているため、着地を待って配る。

---

## 第 1074 波 —— メタ第 12 回を採用した + Y10 が段 2 を閉じた(2026-09-07)

### ★★索引が「偽と判明した statement」を在庫として配っていた（M38。採用済）

メタ第 12 回が **M24 の測定中に別の壊れ方を見つけた**。
★**`decl-index.mjs` は注釈を 1 度も飛ばしていなかった** ——
`.cache/mathlib-index.txt` に **648 行**、`.cache/decl-index.txt` に **38 行**の
「ソースではコメントの中＝実在しない宣言」が載っていた。

★★**これは M14/M20/M22/M24 と向きが逆の害である。** あれらは「在るのに引けない」、
これは**「無いのに引ける」**。ABC3 の 38 件には
`Check/PGC/Theorem42Degenerate.lean:21 theorem theorem_4_2` と
`Check/PGC/FreeTermFunctionRefutation.lean` の `prop_2_2` / `cor_3_1` / `cor_3_3` ——
★**docstring 冒頭に「現行形は偽」と書いてある引用**——が入っていた。
⇒ **索引だけを見た実装者は「`theorem_4_2` は在る」と読む。**

同じ回路で `namespace` / `end` も注釈から拾っていた。実例:
`Topology/AlexandrovDiscrete.lean:40` の docstring の折り返し
「… in the root / namespace instead. -/」を `namespace instead.` と読み、
★**以降 46 宣言に `instead.` が被っていた**。

**本体で検算した(採用後)**:
| 検査 | 結果 |
|---|---|
| `theorem_4_2`(注釈の中)が索引から消えたか | ★**0 件（消えた）** |
| `instead.` の名前空間汚染 | ★**0 件（消えた）** |
| `AlgebraicGeometry.ΓSpec.adjunction` | 0 → **1 件** |
| `AlgebraicTopology.DoldKan.Γ₀.obj` | 0 → **1 件** |
| `lp.evalCLM` | 0 → **1 件** |
| `mathlib-index.txt` | 249,273 → **248,636**(−637) |
| `decl-index.txt` | 24,255 → **24,285**(★新着 3 ファイル分 +67、注釈由来 −38 の差し引き) |
| `--mathlib` 所要 | 4.4 → **5.2 秒**(+0.7) |

★**M22 は既に master に入っていた** —— メタ係は起動直後に「もう直っていないか」を
30 秒で数えて確認し(`._root_.` 3,041 → 0)、**やることが無いと判断して M24 に移った**。
☆**台帳を鵜呑みにせず先に数える**という運用が効いた実例。

★**M24 の実測は台帳の見積り(133 件)の 2.4 倍(325 件)**だった。
台帳は `namespace` 行だけ数えていて **`end` 58 件・`section` 40 件**を数えていなかった。
害は「拾えない」ことではなく**開閉の非対称**(`namespace Isδ₀` は `Is` として積まれ
`end Isδ₀` は落とせない ⇒ 以降の宣言に `Is.` が残る)。

### ★`brief.mjs` の誤報 20 件を直した（M39。本体で実装）

★**`proofParagraphOf` は `kind` が見出し語かを検査していなかった。**
`itemKeyOf`(check.mjs G1 と共有)は `Section` / `Chapter` / `Step` / `Assertion` も拾うので、
`key = "Section 2"` が [pGC] の `.txt` 138 行目の節見出し
「Section 2: Higher Ramification Groups」に当たり、
★**その手前 13 行が 9 件の項目すべての「証明」として出ていた**。

`HEADWORD` 検査 2 行を足した。**退行測定(3 論文 137 項目)**:
Proof が出た件数は **31 / 16 / 4 のまま不変**。
pGC の `noproof` が 29 → 9 に減り(20 件が正しく「原典に見出しを持たない単位」になった)、
Yoshida08 の「見出し不明」17 → 0(同じ 17 件が正しい経路を通るようになった。表示は同じ)。
あわせて `END` に `Q.E.D.` を足した(メタの提案。BK の 6 件。★空文字列に当たらないことを検算済み)。

☆メタ係の測定「`END` に足すべきものは無い」は**正しかった**——
記号が残る 7 本では**見つかった記号の数と `^Proof` 行の数がきっかり一致**する
(Yoshida08 `□` 41/41、GenEll `⃝` 14/14、pGC 3/3、CorrHyp 12/12、EtTh 37/37、
LocProP 55/55、SemiAnbd 36/36)。★**Falt1 は OCR が `□` を小文字 `o` にしているが、
`o` を足してはいけない**("also" "two" "zero" に当たる)。

### ★★本体の見立ての訂正 —— 誤り #7 は「grep で反証できた」は成立しない

本体はメタ係に「7 件目(`pdftotext` は `⊂` を落とす、という誤り)は
**`0_Source/*.txt` を grep すれば 1 秒で反証できた**」と書いた。
★**これは誤りだった。** メタ係の実測:

- pGC の `.txt` の `⊂` は **0 件**なので、grep は**偽の主張を裏づけてしまう**
  (真相は「pGC は `⊂` ではなく `⊆` を使っている」)。
- 反証には **(a) その字を実際に含む対照(Yoshida08)** と
  **(b) `.txt` ではなく `pdftotext` を走らせること**の両方が要る
  (実測: Yoshida08 で pdftotext `⊂` = 53、PyMuPDF = 55)。
- ★**`.txt` は pdftotext 製ではない** —— [[pdftotext-two-implementations-hazard]] の罠の **2 例目**。

⇒ ★**brief の誤り 7 件のうち、機械で捕まえられるのは 2 件、半分が 1 件、4 件は捕まえられない**
(メタ係の結論)。★**「捕まえられない」と正直に測ったことに価値がある。**

### Y10 着地 —— `e(K′/K′′) = |H|` を本当に証明した（Y9 の `hram` が埋まった）

`Found/PGC/FixedRingRamificationIndex.lean` **531 行、13 宣言、sorry 0**。
`#print axioms` 4 本すべて `[propext, Classical.choice, Quot.sound]`。
`lake build ABC3.Found.PGC.FixedRingRamificationIndex` 成功(2987 jobs, 8.3 秒, warning 0)。
往復 **11 回**(うち失敗 5 回、★原因はすべて mathlib の改名)。

★★**本体の段取りは順序が誤っていた。実装者が逆順にして閉じた** ——
本体は「`addVal C ϖ₀ = 1` を示せば `e = |H|` が出る」と書いたが、
★**`addVal C ϖ₀ = 1` を直接示す道は無く、逆順が正しい**:
1. 道 A の前半(`ϖ₀ := ∏_{τ∈H} τ•π` は H 不変・`addVal B ϖ₀ = |H|`)から **`e ≤ |H|`**
2. `M := span_C{π^i : i<e}` に対し「`v(b) ≥ j` なら `M` の元を引いて `v ≥ e` にできる」を
   **j の有限降下帰納**で示し、★**中山の補題**
   (`Submodule.le_of_le_smul_of_le_jacobson_bot`)で `M = ⊤`
   ⇒ `π` は `C` 上**次数 e の monic 関係**を満たす(`exists_pow_eq_sum`、★これが抽象核)
3. 係数が H 不変 ⇒ 相異なる共役 `τ•π` が全部その根 ⇒ 次数比較で **`|H| ≤ e`**
★**`addVal C ϖ₀ = 1` は系として落ちた**(宣言 25 行・証明本体 14 行)。「山」ではなくなった。
★**道 B(Artin + `e·f = n`)は使わなかった** —— Y9 が「mathlib に橋が無い」と測った箇所は迂回できる。

★**抽象核が 16 回連続で効いた**。段 1(`addVal_map_eq_mul_addVal`)は **測定ノイズ以下**、
中山の `exists_pow_eq_sum` が **0.31 秒**、根の数え上げ `card_le_of_pow_eq_sum` が **0.28 秒**。
★**最大の山と目された中山の段は 1 往復で通った。**

★**Y9 の在庫は 4 本中 3 本がそのまま使えた**(`addVal_prod` / `prod_X_sub_C_dvd_of_forall_isRoot` /
`smul_coeff_conjProd` + `eval_prod_X_sub_C`)。
使えなかった 3 本は **Y9 の設定 `SMulCommClass G A B` が Y10 の設定と合わない**ため
(Y10 では `G` は `C` を固定しない。H だけが固定する)。★**設計の差であって欠陥ではない。**

★**本体の指示に無かった退化条件を 1 つ追加発見**: **`[FaithfulSMul G B]` を落とすと偽**
(H が自明作用なら `C = B`, `e = 1`)。
★`e = 0` を許すと段 1 が偽で、しかも**単射性からは出ない**(反例 `ℤ_p ↪ ℚ_p[[t]]`)。
`isUnit_of_isUnit_map` で塞いだ。

★逸脱 4 件。うち 1 番: **`hres`(`K′/K′′` の完全分岐)を結論の側で仮定に置いた** ——
「`K′/K` 完全分岐 ⇒ `K′/K′′` 完全分岐」がこの木にも mathlib にも無いため。
⇒ ☆**新ノード候補**: 上流が Lemma 6.8 を `..._of_fixedRing` に差し替えるなら、
`hres` と `hAC`(`𝒪 ⊆ 𝒪_{K″}`)を供給する小ノードが 1 つ要る。

`lean-idioms.md` **#107(Y10 版)**: 2026-09-07 に名前が動いていた 5 つ
(`Irreducible.not_unit`→`not_isUnit` / `le_or_lt` は無い / `isUnit_of_mul_eq_one`→
`IsUnit.of_mul_eq_one` / `Algebra.algebraMap_mem`→`Subalgebra.algebraMap_mem` /
`IsIntegral` を無名構成子で開くと goal は `aeval` ではなく `eval₂`)。

### Lemma 4.5 着地 —— 抽象核は本当に「1-コサイクル」だった

`Found/PGC/UniformizerCocycle.lean` **519 行、sorry 0**。
`#print axioms` 主要 11 宣言すべて `[propext, Classical.choice, Quot.sound]`
(`IsZCocycle.ext` は `[propext, Quot.sound]` のみ)。
`lake build ABC3.Found.PGC.UniformizerCocycle` 成功(3433 jobs、本モジュール 8.9 秒)。
`lean_start` 1 回(12.3 秒、`lean_status` で imports 一致を確認。差し替わり無し)、
`lean_check` 14 回、`leanfile.mjs` 2 回。

★★**段取りの見立て「抽象核は 1-コサイクル」が当たった。**
`IsZCocycle` / `IsZCocycle.ext`(一意性) / `zpowProd`(存在)が**純群論**で書け、
★**抽象版 Lemma 4.5 本体は `IsZCocycle.ext` に投げるだけで 0.06 秒・6 行**。
★**原典の "argue by induction in both directions" は `IsZCocycle.ext` 1 本に完全に吸収された。**

★**負の `j` が入った**(すべて `j : ℤ`。`ℕ` に逃げていない)。
符号の場合分けは `zpowProd`(`toNat` の一様式)と `zpowProd_step` の証明の 2 箇所に閉じ込めた。
★**第 2 主張は代入で出た**(第 1 主張を `∀ π′, ∀ θ` の一般形で述べたので `π′ := ϕ π`・`θ := π`)。

★**`uniformizerProd` を延ばさず ℤ 版 `zpowProd` を新設した**(測って決めた)。理由:
`uniformizerProd` は `[Monoid G] [MulSemiringAction G B]` の上にあり
(i) `j : ℤ` の `σ^j` が作れない (ii) `B` が環なので `j<0` の `(π_{−j})⁻¹` が作れない。
★`M := Lˣ`・`σ : MulAut M` へ移すと両方が同時に解け、しかも
**分母の非零を `Units.ne_zero` が無条件に供給する**(除算の前提を持ち回らずに済む)。

★★**逸脱 2 は重要**: 「`Θ` の `θ ∈ 𝒪_L`」を落とした。
★**落とさないと原典の第 2 主張が `j < 0` で意味を持たない**(`π_j ∉ 𝒪_L`)。
逸脱 1: `Θ^L` を商でなく**交差積** `θ^ϕ·π = θ·π′` で書いた。
★根拠は**原典 Definition 3.3 自身の "It is an additive group"**
(加法群なら `0 ∈ Θ` で、`0^ϕ/0` は定義できない)。
☆★**この 1 文は `0_Source` の `.txt` の直読で見つかった** ——
`grep -an "Definition 3.3" -A 10`。★**これが無ければ `θ = 0` を除外する不要な仮定を付けていた。**

★**思ったより高かったのは `zpowProd_step`(符号の場合分け)だけ**で、ここだけ 4 往復。
`Finset.prod_range_succ'`(先頭を切る版)が要ることに気づくまでが山だった。
`lean-idioms.md` **#109**(ℤ 上の両方向の帰納法で踏んだ 4 つ。
★**環自己同型の `ℤ` 乗は `RingAut R →* MulAut Rˣ` を 1 本作れば場合分けが消える**)。

### ★字形の凡例について —— 「出なかった」は正しい挙動である

Lemma 4.5 の実装者は「**字形の凡例は私の 2 件の brief には表として出なかった**」と報告した。
★**これは設計どおりである** ——**引用の中に実際に出ている字形だけ**を表にする
(noise を出さないため)。`lemma-4-5` の Proof(401–407 行)には
`Q` も `P` も `̸=` も `b`+大文字も含まれていないので、**出さないのが正しい**。
☆ただし**実装者からは「役に立たなかった」と見える**ので、
**次に凡例が出た持ち場で「役に立ったか」を聞き直す**こと(Prop 4.7 の持ち場に入れた)。

### 新しく配った持ち場: `Proposition 4.7`（Weil 群への延長）

★**Cor 4.9 の後半 `ρ_{f,m} = ρ_{f′,m}` は、その `ρ_{f,m}` 自体が木に無い**ので、
**先に Prop 4.7(ii) が要る**(本体が測って確認した)。

★**在庫が両側から揃っている**のが分かった:
- **ℤ 方向** = 今日 Lemma 4.5 が作った `uniformizerZ` / `zpow_apply_mul_uniformizerZ`
- **単数方向** = 既存の `LubinTateReciprocityIsomorphism.lean`(**425 行**)の
  `galoisUnitReciprocityMap` / `_injective` / `_surjective` / `galoisReciprocityEquiv`
  ——★**これが Yoshida `Prop 4.4(iii)` に当たると本体は見ているが、型で確かめていない**ので
  「違ったらそう報告すること」を持ち場に明記した。

★**未確定の点を持ち場に明示した**: `∐_{j∈ℤ} µ^{(j),×}_{f,m}` の `∐` は、
構造化係が「170dpi では `∐` と `⋃` が紛らわしい」と申告している箇所である
(どちらも `pdftotext` 出力が空)。Lean では `Σ j : ℤ` でも「`j` が `x` から決まる」形でもよく、
**どちらで書いたかを報告**させる。

---

## 第 1075 波 —— ★★Herbrand が形式化された(2026-09-07)

### Y11 着地 —— Proposition 6.9 (Herbrand) + その土台 3 つ、**全部入った**

`Found/PGC/HerbrandFunction.lean` **596 行、宣言 27 本、sorry 0**。
`#print axioms` 主要 7 本とも `[propext, Classical.choice, Quot.sound]`
(`ultrametric_mul_eq_min` は `[propext]` のみ)。
`lake build ABC3.Found.PGC.HerbrandFunction` 成功(2988 jobs, 6.8 秒)。往復 **8 回**。

★★**Λ7 の欠落 3 つのうち 2 つ目(N8 = Herbrand)が閉じた。**
主結論 `herbrand_mem_iff` / `herbrand_coe_mul_coe_eq`(`G_n H = (G/H)_{φ_H(n)}` の集合版)。

★**段取りの見立て「先に立てる節点が 3 つある」が当たり、4 つとも 1 波で入った**:
§1 `i` の超距離不等式(`min_ramIndex_le_ramIndex_mul` / `ramIndex_inv`。★在庫に無かった) /
§2 実数添字の `G_n`(`ramificationGroupReal`) / §3 `φ_H`(`herbrandPhi`) / §4 Herbrand 本体。

★**抽象核が 17 回連続で効いた**。`ultrametric_mul_eq_min` / `exists_max_on_coset` /
`upperSubgroup` / `mem_mul_subgroup_iff` / `RealLeENat` / `truncENat` の 6 本が
★**初回の往復で 0 エラー**(`exact?` 4 本と同居させた 13.58 秒の往復に全部入っていた)。

### ★★原文より安い道が 2 つ見つかった

1. ★**原文の 2 つの場合分け(`i(τ) ≥ m` / `i(τ) < m`)は不要だった。**
   `ultrametric_mul_eq_min` の仮定は `f(στ) ≤ f(σ)`(剰余類の中で σ が最大)の 1 本だけで、
   そこから `f(στ) = min{f(τ), f(σ)}` が 2 行で出る。
   ★原文が場合分けしている理由は `i(σ⁻¹)` を経由する不等式の向きを説明するためで、
   Lean では **`min_eq_right` 1 回で吸収される**。
2. ★**`σ ∈ G_n H/H ⟺ m ≥ n+1` は「最大値が m」だけで出る**(単調性 1 本)。
   代表の取り替えを済ませてしまえば**両向きとも 1 行**。

### ★★除算は `φ_H` の定義式 1 箇所に閉じ込められた（本体の見込みどおり）

★`ℝ` で持った(`ℚ` でも掛け算形でもない)。理由: 原典が `n ∈ ℝ≥0` と**実数の添字**を要求しており、
`(G/H)_{φ_H(n)}` を書くには `φ_H : ℝ → ℝ` でなければならない。
`ℚ` にすると `G_n` の添字型と食い違い、`ℚ ↪ ℝ` の橋が余計に要る。

★★**証明本体は `card_mul_herbrandPhi_add_one`(`|H|·(φ_H(n)+1) = herbrandSum`)だけを使い、
除算を一度も見ない。** ⇒ ★**Y9 の掛け算形 Lemma 6.8 がそのまま噛み合った。**
☆**#102(「`ℕ∞` の除算・切り詰め引き算を書かない」)以来の方針が、
「除算が原理的に避けられない場所」でも定義の 1 行に閉じ込める形で生き延びた。**

`min{i(τ), n+1}` は `ℝ≥0∞`/`EReal` を使わず自前の
`truncENat (x : ℕ∞) (r : ℝ) : ℝ := if x = ⊤ then r else min (x.toNat : ℝ) r` に落とした
(理由: 2 段 coercion より `cases x` + `simp` の方が補題名を探さずに済む。関連 5 本はすべて 1–3 行)。

### ★添字の一致を確かめた（下流がずれないため）

`realLeENat_natCast_iff : RealLeENat (m:ℝ) x ↔ (m:ℕ∞) ≤ x` を橋にして、
`ramificationGroupReal_eq_of_mem_Ioc`(★**原典の "`G_n = G_i` if `i ∈ ℤ≥0` and `n ∈ (i−1,i]`"
そのもの**)と系 `ramificationGroupReal_natCast` を証明した。
★**原典が `(i−1, i]` と書いている両側の境界を両方仮定に置いてあるので、
`n ≥ 0` を落としても下流がずれない**(`n > −1` が自動で出る)。

### ★Y10 の成果がそのまま使えた

`card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing`(Y10 の `hram` 無し版)は
★**引数順を写すだけで通り、本ファイルが新たに足した仮定は 1 つも無い**。
`hres`(`K′/K′′` 完全分岐)と `hAC` は依然として仮定のまま。

### ★★`leanfile.mjs` の測定限界が判明した（D27 の逃げ道の唯一の実質的な損）

★Y11 の実測: `leanfile.mjs` は 1 往復 **9.6–13.6 秒**で、
**うち import 読み込みが 10.7 秒**(`example : 1 = 1 := rfl` だけのファイルで実測)。
⇒ ★★**宣言ごとの秒数(抽象核 0.05 秒 / 具体層 2 秒 の刻み)は測定ノイズに埋もれる。**
**MCP `lean_check` の 0.05 秒刻みは leanfile では再現できない。**
☆これが「MCP は 1 体・2 体目は leanfile」という D27 の規約の**唯一の実質的な損**である。
⇒ ★**leanfile 側の持ち場では「往復何回・一発かどうか」を報告させることにした**
(秒数を求めると嘘を書かせることになる)。

### ★brief の 1b は 5 人目も「決定的」と報告

Y11:「1b の 6 行がそのまま §1〜§4 の設計図になり、`lean-search` を 1 回も呼ばずに済んだ。
★**字形の凡例(`φ` が U+03C6、`Σ` と分数の横線が pdftotext に出ない)のおかげで
docstring の逐語を作り直す往復が 0 回**だった」。
⇒ ☆**凡例は「出た持ち場」では効いている**(Lemma 4.5 では該当字形が無く出なかった)。

`lean-idioms.md` **#110**(`Subgroup G` を `Set G` に落とすとき `G` が推論できない)・
**#111**(`field_simp [f]` は `def f` を展開しない)・
**#112**(★`push_neg` が 2026-09-07 に deprecated)。

### 新しく配った持ち場: `Lemma 6.10`（Hasse-Arf の直前）

★本体が原文(`.txt` 1077–1101 行)を直読し、**(ii) が独立の山**であることを確認した:
「φ は連続で区分線形だから**微分だけ比べればよい**」という 3 行が、Lean では
**「区分線形関数の一致判定」という節点**になる(構造化係も同じ指摘をしていた)。
⇒ ★**2 つの道(A: 原文どおり微分 / B: 差分。微分積分を使わない)を挙げ、
どちらが安いかを測って選ばせる**ことにした。★本体は確かめていない。
★Y11 の `ramificationGroupReal_eq_of_mem_Ioc` が**区間ごとの議論の土台として既に在る**。

★**Hasse-Arf (Thm 6.11) にはさらに 2 つ要る**(本体が原文 1102–1124 行を直読して確認):
**Prop 6.2**(`G_1` が p 群・`G ≅ ⊕ ℤ/p^{m_i}ℤ`)と **Cor 6.3**(`G_n ≠ G_{n+1} ⇒ e_0 ∣ n`)。
★**どちらもまだ木に無い見込み**(★Cor 6.3 は索引で未確認)。

### ★Hasse-Arf に要るものを索引で測り直した —— **欠けているのは Cor 6.3 だけ**

本体は前波で「**Prop 6.2 と Cor 6.3 はどちらもまだ木に無い見込み**(★Cor 6.3 は索引で未確認)」と
書いた。★**測ったら半分が誤りだった。**

| 原典 | 木の状態 | 実体 |
|---|---|---|
| **Prop 6.2**(θ_0・θ_n の 2 本の単射) | ★**在る** | Y3 `RamificationQuotientEmbedding.lean`(745 行)の `ramCoeff` 一族 |
| `G_1` が p 群 | ★**在る** | `isPGroup_lowerRamificationGroup_one` / `_Adjoin` |
| **Cor 6.3**(`G` 可換 ∧ `G_n ≠ G_{n+1}` ⇒ `e_0 ∣ n`) | ★★**無い** | —— |
| Cor 6.7 | 在る | Y8 `SenJumpExpansion.lean` |
| Prop 6.9 (Herbrand) | 在る | Y11 `HerbrandFunction.lean` |
| Lemma 6.10 | **走行中** | Y12 |
| 有限アーベル p 群の構造定理 | 未測定 | mathlib にあるはず(★型で引くこと) |

⇒ ★**Hasse-Arf の手前で本当に欠けているのは Cor 6.3 の 1 本だけ**である。

### ★Cor 6.3 は在庫が揃っている（本体が原文と型の両方を読んだ）

原文(`.txt` 930–945 行、`□` まで)の骨格:
1. `θ_n(στσ⁻¹)` を `π′ = σ⁻¹(π)` で計算する。`τ(π′) = π′(1+a)`(`a ∈ 𝔭^n`)なら
   `στσ⁻¹(π) = π(1 + σ(a))` ゆえ `θ_n(στσ⁻¹) = σ(a) mod 𝔭^{n+1}`。
2. `a = bπ^n`・`σ(π) = uπ` と書くと `σ(a) = σ(b)u^nπ^n ≡ b u^n π^n = u^n a (mod 𝔭^{n+1})`
   (`σ(b) ≡ b mod 𝔭`)。⇒ ★**`θ_n(στσ⁻¹) = u^n · θ_n(τ)`**。
3. `G` 可換なら `στσ⁻¹ = τ` ゆえ `a ≡ u^n a`。`G_n ≠ G_{n+1}` から `θ_n(τ) ≠ 0` なる τ が取れ、
   `θ_0(σ) = u` が `G_0/G_1` を生成するよう σ を取れば **`e_0 ∣ n`**。

★**在庫(型で引いて実測)**:
- ★**`ramCoeff_smul_conj : ramCoeff (τ • α) n (τ*σ*τ⁻¹) = τ • ramCoeff α n σ`**
  (`RamificationJumpDivisibility.lean:164`) —— ★**上の 1.+2. がほぼそのまま在る**
- `ramCoeff_independent`(`RamificationQuotientEmbedding.lean:480`) —— π 非依存性
- `ramCoeff` 一族 10 本(`ramCoeff_mul` / `pow_mul_ramCoeff` / `isUnit_one_add_ramCoeff_zero` ほか)

⇒ 見積 **250–500 行**。★**実装枠(2 体)が空き次第これを配る。**

★★**記法の罠**: 原典の `G_n ≠ G_{n+1}` は `pdftotext` で **`=` にしか出ない**
(斜線がベクター描画)。★**テキストだけを読むと主張が反転する**
([[pdftotext-drops-negation]])。構造化係は `data-txt="="` で明示済み。
★**持ち場に必ず書くこと。**

### Prop 4.7 着地（部分）—— (i) の全単射と (ii) の同型性が入り、4 件を新ノードに切り出した

`Found/PGC/WeilReciprocityExtension.lean` **697 行、宣言 42 件、sorry 0**。
`#print axioms` 12 宣言すべて `[propext, Classical.choice, Quot.sound]`
(`negFiberSigmaEquiv` と `weilDegKerEquiv` は **`Classical.choice` すら不要**)。
`lake build ABC3.Found.PGC.WeilReciprocityExtension` 成功(3434 jobs、当該モジュール 7.7 秒)。
MCP 側なので秒数が測れた: 抽象核 **0.25–0.55 秒** / 具体層 **0.26–0.34 秒**、
`lean_check` の総計 **3 秒程度**。

| | 入った | 落とした |
|---|---|---|
| (i) 全単射 | ★`prop_4_7_i_equiv` | — |
| (i)「`L^m_f` は `K` 上 Galois」 | — | ★**全部**(Prop 4.4(ii) と Lemma 4.6 が木に無い) |
| (ii) `ρ` が同型 | ★`weilReciprocityEquiv` | — |
| (ii) `ρ` の構成 | 準同型性の代数計算のみ | ★`ρ` 自体(`µ_{f,m}` 上の `W` の作用が要る) |
| 極限 `ρ_f` | — | ★全部 |

### ★★`∐` の読みが確定した —— 構造化係の申告が解決した

構造化係は「170dpi では `∐` と `⋃` が紛らわしい」と申告していた(どちらも pdftotext 出力が空)。
★**実装者が Lean で確かめて `∐` が正しいと確定させた**:
`Σ j : ℤ, B j`(依存和)で書くと型として非交和になり、
`negFiberSigmaEquiv_fst : (Θ x).1 = -w x` が「`j` は `x` から `v(x) = −j` で決まる」を与える。
★★**`⋃` と読むと `j` の一意性が落ちて像の合併になり、全単射が主張できない。**
⇒ 原典の丸括弧 `(v(x) = −j)` が `j` を一意に決めているので **`∐` が正しい**。

### ★★本体の推測が誤っていた（今日 8 件目）—— ただし「未確認」と書いたので覆った

本体は持ち場に「`galoisUnitReciprocityMap` 系(既存 425 行)は Yoshida `Prop 4.4(iii)` に
当たると見ているが、**型で確かめていない**ので違ったら報告すること」と書いた。
★**実測の結果、違った。**

`galoisReciprocityEquiv` の型は
`(K.carrier⟮x⟯ ≃ₐ[K.carrier] K.carrier⟮x⟯) ≃* (𝒪[K.carrier] ⧸ Ideal.span {π^n})ˣ`
すなわち ★**`L = K`(p 進局所体そのもの)の場合の古典的 Lubin-Tate 主定理**である。
Yoshida の Prop 4.4(iii) は `L` = `K` の**完備不分岐拡大**に対するもので、
★★**Prop 4.7(ii) が使うのは `L = K̂`(不分岐閉包の完備化)の場合**。
⇒ ★**既存 425 行は Prop 4.7(ii) の入力を直接は供給しない。**
さらに `L^m_f = L(µ_{f,m})` ではなく `K.carrier⟮x⟯`(原始点 1 個の添加)で書かれており、
その同一視も statement に無い。

★実装者は `u : w.ker ≃* Rˣ` / `e : res.ker ≃* v.ker` を**仮定に取る**ことで、
`K̂` 版が立った時点で代入できる形にした。★**通すために statement をねじ曲げていない。**

☆★**「推測を推測として渡す」が 2 度目も機能した**(1 度目は #8 SubfieldClosed の 45 行)。

### ★★字形の凡例が決定的だった（2 件目の証言）—— しかも設計を救った

Prop 4.7 の実装者:
> 字形の凡例は今回は表として出た(`b + 大文字` → ハット、`∼ =` の 2 行割れ → `≅`)。
> ★**これが決定的で、`.txt` の `L = bK` が `L = K̂`(K ではない)だと分かった**
> —— これが無いと「(ii) は `L = K`」と読んでしまい、`W(K̂^m_f/K)` の設計を丸ごと間違えていた。

★**Lemma 4.5 で凡例が出なかったのは「引用に該当する字形が無ければ出ない」設計どおりで、
正しい挙動だと確認できた**(実装者自身がそう書いた)。
⇒ ☆**前波で積んだ「次に凡例が出た持ち場で聞き直す」が解決した。**

### ★Weil 群は「対」のまま持った（本体の見立てどおり安かった）

`weilGroup res ϕ : Subgroup (Γ × Multiplicative ℤ)`、carrier は `{(σ,j) | res σ = ϕ^j}`。
理由: (a) 群構造が積の部分群として自動で付く(**0.25 秒**) (b) `deg` が第 2 成分の射影なので
`Classical.choose` が不要 ⇒ ★**`ϕ` の位数無限を仮定しなくてよい**
(`Subgroup.zpowers ϕ` の逆写像を作る道だと無限位数が必須になる)
(c) 原典の元の書き方(`ϕ^j on K̂, α ↦ [xπ_j](α)`)とそのまま一致する。

★**Lemma 4.5 の投資がそのまま効いた**: ℤ の両方向帰納法(#109)は `Int.induction_on` の
3 場合を 1 回書いただけで済んだ(`uniformizerZ` 側が全部 `j : ℤ` で用意されていたため)。

`lean-idioms.md` **#113**(ℤ で添字づけた `Σ` 型を作るときの 4 つの穴。
★`Equiv.sigmaFiberEquiv` と `Equiv.sigmaCongrLeft'` の合成は `Eq.mpr` が挟まって
`_apply` が `rfl` にならない、ほか)。

### 新しく必要になったノード（Prop 4.7 が切り出した 4 件）

1. **Prop 4.7(i) 前半「`L^m_f` は `K` 上 Galois」**(Prop 4.4(ii) と Lemma 4.6 が前提)
2. **Lemma 4.6**(`[θ] : µ_{f,m} ≅ µ_{f′,m}` かつ `L^m_f = L^m_{f′}`)
3. ★**Proposition 4.4 の (ii)(iii) を `L = K̂` で**(既存 425 行は `L = K` 版なので流用不可)
4. **`ρ_f`(`m → ∞` の極限)**

### 新しく配った持ち場: Y13 `Corollary 6.3`（Hasse-Arf の最後の欠落）

★在庫が揃っていることを本体が型で確かめてから配った(`ramCoeff_smul_conj` /
`ramCoeff_independent` / `ramCoeff` 一族 10 本)。
★★**`ramCoeff_smul_conj` が原文の `θ_n(στσ⁻¹) = u^n θ_n(τ)` に当たるという見立ては
実機で確かめていない**ので、「違ったら遠慮なく訂正して報告すること」を明記した
(★本体の誤りが今日 8 件目になったことも書いた)。
★**記法の罠を最重要事項として書いた**: 原典の `G_n ≠ G_{n+1}` は
`pdftotext` に `=` としか出ず、**テキストだけ読むと主張が反転する**。

### Y12 着地 —— Lemma 6.10 (i)(ii) 両方、★★第 3 の道で閉じた

`Found/PGC/HerbrandComposition.lean` **699 行 / 29 宣言、sorry 0**。
`#print axioms` 10 宣言すべて `[propext, Classical.choice, Quot.sound]`。
`lake build ABC3.Found.PGC.HerbrandComposition` 成功(2989 jobs, 7.1 秒)。
`leanfile.mjs` 往復 **13 回**、MCP は指示どおり一度も触っていない。

★**落とした主張は無い。** 逆に (i) は原典の 2 文(`φ_G(0)=0` と `n ∈ ℤ≥1` の閉じた形)が
`n : ℕ` の**1 本の式に収まった**(`n=0` で `Finset.Icc 1 0 = ∅`)。★**原典より強い。**

### ★★(ii) は道 A（微分）でも道 B（差分）でもなかった

本体は 2 つの道を挙げたが、実装者は**どちらでもない道**で閉じた:
掛け算形で書くと (ii) は `Σ_{σ∈G} min{i(σ),n+1} = Σ_{σ∈G} min{i_ϖ(σ), φ_H(n)+1}` と同値で、
★**`H` の剰余類ごとに `min` の結合則だけで一致する**。
剰余類 `σH` の最大代表 `σρ`(`m := i(σρ)`)について両辺とも `Σ_{τ∈H} min{i(τ), min{m, n+1}}` になる。

- ★**`deriv` も `Continuous` も使っていない。** 解析的な道具は `Monotone.map_min` **1 本だけ**。
- ★★**「mathlib に区分線形関数の一致判定が在るか」は調べていない ——
  探す前に不要になった。** ☆**実装者はこれを正直に書いた**(在庫照会をサボったのではなく道が消えた)。
- ★★**剰余類への分割すら要らなかった**: `Σ_{σ∈G} Σ_{τ∈H}` を `Finset.sum_comm` で入れ替え、
  `σ ↦ στ` が `G` の全単射(`Equiv.mulRight`)であることを使うと `|H| · Σ_{σ∈G}` になる。
  ⇒ ★**商群 `G ⧸ H` を作らずに閉じる**(Y9・Y10・Y11 と同じ流儀が 4 波続いた)。

### ★`|G_n H/H|` を取り出す必要が生じなかった

原文の `G_n/H_n = G_n/(H ∩ G_n) ≅ G_nH/H`(**第 2 同型定理**)は**この道では一度も現れない**。
★難しくて避けたのではなく、微分を比べないので `|(G/H)_{φ_H(n)}|` も `|G_nH/H|` も出てこない。
⇒ ★**Y11 が挙げていた「`(G_n : Set G) * (H : Set G)` を `Subgroup` に述べ直すノード」は
不要になった**(必要なら別途立てればよいが Lemma 6.10 のためには要らない)。

代わりに `φ_{G/H}` の同定を `phiOf_comp_of_card_fiber` → `phiOf_quotient` →
`herbrandPhiGroup_eq_phiOf_quotient` で立てた。根拠は `|G| = |H|·|G/H|` と
「`i_ϖ` は `H` 剰余類上で定数」の 2 つだけ。

### ★★添字ずれ（`Σ_{i=0}` 対 `Σ_{i=1}`）を構造的に防いだ

★**抽象核を 2 段に分けた**: `phiOf_natCast` が `Σ_{i=0}^{n}`(**原文の途中式そのもの**)、
`phiOf_natCast_of_pos` が `f > 0` を使って `i=0` の項(`= |S|`)を `−1` と相殺して
`Σ_{i=1}^{n}`(原文の結論)。
★★**2 つの差がちょうど `pos_ramIndex`(`i(σ) ≥ 1`)であることを退化検査に書いた。**
☆この論文には既に erratum が 2 件あるので、**ずれる場所を宣言の形で見えるようにしたのは良い設計**。

### ★配管の実測 2 つ

- ★**`Fintype ↥(⊤ : Subgroup G)` は `[Fintype G]` から推論されない**(probe で確認)。
  ⇒ `φ_G` は `G` 全体の上の和として `herbrandPhiGroup` を新規定義し、
  `herbrandPhiGroup_eq_herbrandPhi_top` で Y11 の `herbrandPhi α ⊤` と一致することを確かめた。
- `lean-idioms.md` **#114**: ★**束縛子の型を書かないと `(τ : G)` が coe ではなく型注釈に読まれる**。
  `fun τ => pos_ramIndex huni (τ : G)` と書くと `τ : G` に固定されて
  **`failed to synthesize Fintype G`** が出る(インスタンスの問題に見えるが原因は elaboration の順序)。
  直しは `fun τ : H => …`。

★**逸脱 3 が重要**: **(ii) で `H` の正規性 `G ▷ H` を仮定していない**。
`φ_{G/H}` を `G` 上の和で書くので商群を作らず、計算は正規性を使わない
(`[H.Normal]` は `φ_{G/H}` の同定補題 3 本にだけ置いた)。

★**見積との差**: 本体は「(ii) は重いから (i) だけでもよい」と書いたが、
★**(ii) の本体(`herbrandPhiGroup_comp`)は 1 往復で一発**だった。
高かったのは**(i) の `ℕ∞ ↔ ℝ` の橋**の方(cast まわりで 4 往復)。
★brief の 1b は「特に (i) の `Σ_{i=0}^{n}|G_i|` という**途中式**が原文に明示されていたので、
添字ずれの罠を最初から避けられた」。字形の凡例(`∼ =` → `≅`、合成記号 `◦` = U+25E6)も効いた。

★★**Hasse-Arf (Thm 6.11) は `herbrandPhiGroup_natCast` と `herbrandPhiGroup_comp` の両方を
そのまま入力にできる形になっている。**

### 新しく配った持ち場: Y14 `hres` / `hAC` の供給（3 ノードが引きずる債務）

★**Y10・Y11・Y12 の 3 本(計 1,826 行)が同じ 2 つの仮定を持ち回っている。**
Y10 と Y11 の実装者が**独立に「供給する小ノードが 1 つ要る」と挙げた**ので配った。
★本体の見立て「`A → C → B` の塔が立てば `hres` はほぼ自明で、実質は塔を立てること」は
**実機で確かめていない**ので、「外れたらそう報告すること」を明記した。
★`C` の構成は道 A(`FixedPoints.subring`)/ 道 B(抽象のまま存在を述べる)を挙げて測らせる。
★**既存 3 本から `hres` を消して回らないこと**を明示した(差し替えはゲート後の本体判断)。

### ★★★本体の手順の穴 —— **道具が出していた警告を自分で切り落としていた**

Y13 の持ち場は「Cor 6.3 は木に無い」という前提で書いたが、★**既に埋まっていた**
(`RamificationJumpDivisibility.lean`、第 1071・コミット `08e1bdbe`、`Found.lean:1745` で import 済み。
`nat_card_quot_zero_dvd_of_ne` / `nat_card_quot_top_dvd_of_ne` /
`nat_card_quot_lowerRamificationGroupAdjoin_one_dvd` の 3 形)。

★★**`brief.mjs` は冒頭 8 行目でちゃんと警告していた**:
```
★**この項目は既に木にある**: `Found/PGC/RamificationJumpDivisibility.lean`
  ⇒ 新規ファイルではない。`node tools/brief.mjs --node …` の方が情報が多い。
```
★★**本体が `sed -n '/## 1\./,/…/p'` で出力を絞ったため、この行を切り落としていた。**
(索引の grep でも `e_?0` 等の語で引いたが、実際の宣言名は `nat_card_quot_*` だったので当たらなかった
——★**「在庫は型で引く」を本体自身が守っていなかった**。)

⇒ ★★**規則: `brief.mjs --paper` の出力は必ず冒頭から読む。`sed` で絞るなら `head` を併用する。**
☆**道具が正しく警告を出していたのに人間側の絞り込みで消していた**、という失敗形である。
★次の持ち場(Thm 6.11)では冒頭を確認してから配った(木に無いことを確かめた)。

★**同じ持ち場でもう 1 つ誤っていた**: 本体は
「`ramCoeff_smul_conj` が原文の `θ_n(στσ⁻¹) = u^n θ_n(τ)` に当たる」と書いたが、
★**違った**。`ramCoeff_smul_conj` は**素元取り替えと共役の整合だけ**で `u^n` の因子を持たない。
原文に当たるのは同ファイル **185 行の `residue_ramCoeff_conj`**:
`residue (ramCoeff α n (τσ6τ⁻¹)) = residue (1 + ramCoeff α 0 τ)^n * residue (ramCoeff α n σ)`。
★**捻れ `u^n` は `residue_ramCoeff_change_uniformizer` 側から来ており、
`ramCoeff_smul_conj` からは来ない。本体の対応表は 1 リンクずれていた。**
☆ただし**持ち場に「実機で確かめていない。違ったら訂正して報告すること」と書いたので覆った**
(この形式が機能したのは 3 度目)。

### Y13 着地 —— 持ち場を組み替えて別の価値を出した

`Found/PGC/AbelianJumpDivisibility.lean` **487 行 / 21 宣言(証明本体は約 90 行)、sorry 0**。
`#print axioms` 9 本すべて `[propext, Classical.choice, Quot.sound]`。
`lake build` 成功(2981 jobs, 6.6 秒)。★**1 度も詰まらず `lean_check` 往復 11 回で全部通った。**

★**既に埋まっていると分かった後、実装者は中身を 2 つに組み替えた**:
- **(A) 原文の段取りに忠実な別証明**(先行ノードは `Monoid.exponent` 経由で原文の経路ではない)。
  `e_0` の 2 つの表し方(`Nat.card (G ⧸ G_1)` と「生成元 `σ` の `θ_0(σ)` の位数」)の
  **一致を証明**した。
- ★★**(B) Hasse-Arf 最終段が Cor 6.3 から取り出す形**(原文が **`hence` の 1 語で畳んだ部分**)。
  `dvd_sum_natCard_lowerRamificationGroup` と `dvd_of_sum_eq_natCard_mul`。
  ★**原文の "hence `e_0 ∣ Σ|H_i|`" は純算術だった**と見抜き、
  `sum_Ico_block_eq` / `dvd_sum_Ico_mul_of_step` / `dvd_sum_Icc_of_step` に切り出した
  (**0.01–0.15 秒**)。

★**`G` 可換を 1 行に閉じ込めた**(`hconj : τ * σ * τ⁻¹ = σ` ただ 1 行)。
★**`G` 可換を `∀ x y : G, x * y = y * x` の命題で渡す**(`CommGroup` 構造だと具体層に当たらない)。
★**逸脱 4 が良い設計**: 最後の一歩(`e_0 ∣ φ_H(n)`)は Lemma 6.10(i) を
`hk : Σ = |G_1| * k` という**仮定として受け取る** ——
Λ7 が Lemma 6.10 を**同時走行中**だったので**依存を作らず「差し込むだけ」の形**にした。

`lean-idioms.md` **#115**(4 項目)。とくに (b) ★**`Nat.mul_le_mul_left e (Nat.le_succ m)` は
`e * m.succ` を作り `omega` が `e * (m+1)` と別原子として扱う**(#103 と同型)。

### 新しく配った持ち場: Y15 `Theorem 6.11 (Hasse-Arf)` —— ★★入力が全部揃った

★**本体は今回、`brief.mjs` の冒頭を確認してから配った**(木に無いことを確認済み)。

| 原典 | 木の宣言 |
|---|---|
| Cor 6.3 | `nat_card_quot_*`(Y4) + ★**Y13 の受け渡し口 `dvd_of_sum_eq_natCard_mul`** |
| Cor 6.7 | `exists_padicDigits_card_lowerRamificationGroup`(Y8) |
| Prop 6.9 | `herbrand_mem_iff`(Y11) |
| Lemma 6.10 (i)(ii) | `herbrandPhiGroup_natCast` / `herbrandPhiGroup_comp`(Y12) |
| Prop 6.2 | `residue_ramCoeff_conj` ほか(Y3) |
| `G_1` が p 群 | `isPGroup_lowerRamificationGroup_one` |

★**未確認として残したのは 1 点**: 有限アーベル p 群の構造定理が mathlib に在るか
(★**本体は名前を確かめていない**)。持ち場には
「無ければ『非自明な巡回商 `G/H ≅ ℤ/p^mℤ` が取れる』だけで回すこと ——
★**原文が実際に使っているのはそれだけ**」と書いた。

### Y14 着地 —— `hres` / `hAC` が本当に埋まった（§3 まで入った）

`Found/PGC/FixedRingTower.lean` **434 行、sorry 0**。
`lake build ABC3.Found.PGC.FixedRingTower` 成功(2990 jobs、自ファイル 6.5 秒)。
`lean_check` 往復 **11 回**。
★`#print axioms` で **`exists_sub_mem_of_exists_image` は "does not depend on any axioms"**。

★★**本体の見立て「塔が立てば `hres` はほぼ自明」は当たった。しかも塔すら要らなかった。**
`hres` の中身は `hAC` から **3 行**(`exists_sub_mem_of_exists_image`)で、
`Algebra A C` のインスタンスを立てる必要がない。
★**本ノードの実質は `hAC` の側**で、それは `[SMulCommClass G A B]` から
```lean
smul_algebraMap_eq_self : ρ • algebraMap A B a = algebraMap A B a := by
  rw [Algebra.algebraMap_eq_smul_one, smul_comm, smul_one]
```
の **1 行**だった。⇒ ★**§1・§2 は本体の見立てより安く、代わりに §3 が重心になった**(約 100 行)。

### ★★`C` を `FixedPoints.subring` で作った理由が実測で出た（道 A）

`fixedRing B H := FixedPoints.subring B ↥H`。★**`Subring B` なので
`Algebra ↥(fixedRing B H) B` / `IsScalarTower` / `IsDomain` が
すべて mathlib のインスタンスで即座に付いた**(`#synth` **0.03 秒**で確認)。
★**道 B(抽象のまま存在を述べる)を採ると `[IsDiscreteValuationRing C]` を
永久に仮定のまま持ち回ることになり、Y10 の債務が消えない。**
☆Y9・Y10 が一貫して採ってきた「抽象のまま受ける」流儀を、**ここで初めて逆に振った**
——**債務を閉じるノードでは具体的に作るほうが正しい**、という設計判断である。

### ★★§3 の発見: `B^H` が DVR であるのに**完全分岐は要らない**

`instance fixedRing_isDiscreteValuationRing`(仮定は `[Fintype ↥H]` のみ)。
★**`H` が有限でありさえすれば `B^H` は DVR** である。
抽象核は「離散付値環 `B` の部分環 `C` が (i) 単元が降りる (ii) 割り算で閉じている
(iii) 非零非単元を持つ、なら `C` は DVR」で、★**群も Galois も分岐も出ない**
(`v(C∖0) ⊆ ℕ` が差で閉じた部分モノイド → `Nat.find` で最小正元 → 強帰納法 →
`ofHasUnitMulPowIrreducibleFactorization`)。**0.47 秒・2 往復**(誤りは `Associated` の向き 1 箇所)。

**主結論**(仮定は**底の設定だけ**):
```lean
addVal_map_eq_card_mul_fixedRing (hπ) (hresA) (hadjA) (c) :
  addVal B (algebraMap (fixedRing B H) B c) = (Nat.card H : ℕ∞) * addVal (fixedRing B H) c
```
系 `addVal_map_uniformizer_eq_card_fixedRing`(`e(K′/K″) = |H|`)。

### ★同時実行での良い判断（記録に値する）

Y14 は `lean_status` の imports が自分の要求と一致しなかったが、
★**自分が引く在庫が全部その部分集合だと確かめて `lean_start` を 1 回も呼ばなかった**
(他体が使っている env を壊さないため)。最終確認だけ `leanfile.mjs` + `lake build`。
☆**D27 の規約が想定していなかった第 3 の選択肢**であり、**次の持ち場に前例として書いた**。

`lean-idioms.md` **#116**: ★(a)**`def` が包む部分構造は担い手の型を明示引数にしないと
メタ変数になり、症状が `MulSemiringAction G ?m` / `OfNat (↥(fixedRing H)) 0 is stuck` として
5 箇所に散る**(原因は型クラスではなく暗黙引数の設計)/ (b) `associated_one_iff_isUnit` の向き /
(c) `rw [← h]` は等式ゴールの左右両方を書き換える。

### 新しく配った持ち場: Y16「`H ⊴ G` なら `G` は `B^H` に環作用する」

★Y14 が供給**できなかった**のは 3 つだけ:
`[MulSemiringAction G C]` / `hcomp` / `hHtriv`。理由は
★**`G` が `B^H` に作用するのは `H ⊴ G` のときだけ**で、Y14 は正規性を仮定していなかったため。
★**実装者自身が「必要なら 1 ノードを別に立てるのが正しい切り方」と提案した**ので配った。

★★**これが入ると Λ7 の鎖(Y10 → Y11 Herbrand → Y12 Lemma 6.10 → Y15 Hasse-Arf、
計 2,300 行超)が原典 §6.1–§6.2 の設定だけに載る。**
★mathlib に既にある可能性があるので「**在ったら在ったと報告すること**」を明記した。

### Y16 着地 —— Λ7 の債務が 8 つ中 7 つ消えた（213 行、思ったより安い）

`Found/PGC/FixedRingAction.lean` **213 行、sorry 0**（見積 150–350 の内）。
`lake build` 成功(2991 jobs, 6.1 秒)。
★`lean_check` **4 回**(0.01–0.61 秒)、`leanfile.mjs` **1 往復で ok(警告 0)**。

★**抽象核が 22 回連続で効いた。しかも段取りより 1 段一般化された**:
```lean
smul_mem_fixedPoints_of_normal {α G} [Group G] [MulAction G α] {H : Subgroup G}
  (hH : H.Normal) (ρ : G) (ha : ∀ τ ∈ H, τ • a = a) (hτ : τ ∈ H) : τ • (ρ • a) = ρ • a
```
★**環すら出てこない**（`MulAction G α` だけ）。**0.14 秒・一発**。
`#print axioms` は **`[propext]` のみ**。

★**段取りの見立てがすべて当たった**（この波は珍しく訂正が無い）:
- `hcomp` は **`rfl` 1 語**、`hHtriv` は **`Subtype.ext` 1 語**（見立て「1–3 行」どおり）
- `H ⊴ G` を使ったのは **1 箇所だけ**（`hH.conj_mem'`。grep で確認済み）
- `[Fintype ↥H]` は §1・§2 では **要らなかった**（§3 でだけ要る）

### ★★mathlib には「商群版しか無かった」

`Mathlib/RingTheory/Invariant/Basic.lean:98` に
`instance (H : Subgroup G) [H.Normal] : MulSemiringAction (G ⧸ H) (FixedPoints.subring B H)`。
★**Y10 が要求するのは `MulSemiringAction G C`(`G` そのもの)なのでそのままでは使えない。**
加えて `Mathlib.RingTheory.Invariant.Basic` は本木のどこからも import されていない
(`Algebra.IsInvariant` が `Unknown constant`、0.01 秒で確認)。
`MulSemiringAction.compHom (QuotientGroup.mk' H)` で引き戻す道はあるが
**mathlib import を 1 本増やすだけで `hcomp` はどちらでも `rfl`** なので直接構成した。
☆**「在ったら在ったと報告すること」と書いた結果、"部分的に在った" という正確な答えが返った。**

★**主結果 `card_mul_ramIndex_eq_sum_ramIndex_fixedRing`**:
Y10 が持ち回っていた `[MulSemiringAction G C]` / `hcomp` / `hHtriv` / `hinj` / `hfixC` /
`hres` / `hAC` / `[IsDiscreteValuationRing C]` が**すべて消えた**。

`lean-idioms.md` **#117**:
★(i)**素の `def` 包み(`fixedRing B H := FixedPoints.subring B ↥H`)は `instances` 透明度で
展開されないので mathlib のインスタンスが降りてこない**(害)が、逆に
**自前インスタンスが mathlib 側と衝突しない**(益)/
★(ii)**在庫調査は名前空間で `.cache/mathlib-index.txt` を 1 回 grep するのが最速**
(`FixedPoints\.` で 60 行、「商群版しか無い」が確定)/
★(iv)★★**この環境では `cat > f << 'EOF'` のヒアドキュメントが PreToolUse フックで壊れる。
Lean ファイルは Write ツールで書くこと**(実際に踏んだ)。

### 新しく配った持ち場: Y17 `hfix`（Λ7 に残った最後の 1 つ）

残った唯一の仮定は
`hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A {ι ϖ}`
= ★**原典 Lemma 5.11 を `K″/K` に当てた形**(`𝒪_{K″} = 𝒪[π″]`)。
★**数学が足りない側の穴であり、配管の穴ではない**(Y16 の実装者の判定)。

★本体が**型で**在庫を測ってから配った(`adjoin_uniformizer_eq_top` /
`exists_algebraMap_fixedRing` / `exists_sub_mem_fixedRing` / `adjoin_fixedRing_eq_top` /
`adjoin_eq_top_of_image`)。
★**本体が特定した重心**: `Algebra A C` のインスタンスが無い可能性が高い
(Y14 は「立てる必要がない」と書いて**避けている**)。★**実機で確かめていないので測らせる。**
★★**`brief.mjs --id lemma-5-11` の冒頭を確認するよう明記した**
——**本体は今回それを確認していない**ので、agent に確認させる形にした
(2026-09-07 の失敗形 [[brief-header-carries-the-warning]] の再発防止)。

★これが埋まると **Λ7 の鎖(計 2,500 行超)が原典 §6.1–§6.2 の設定だけに載る。**

### Y15 着地（部分）—— Hasse-Arf の 3 段のうち **2 段が sorry 0 で埋まった**

`Found/PGC/HasseArf.lean` **615 行、sorry 0**。
`#print axioms` 全部 `[propext, Classical.choice, Quot.sound]`。
`lake build ABC3.Found.PGC.HasseArf` 成功(2991 jobs, 6.8 秒)。`lean_check` **13 往復**(うち失敗 3)。
`lean_start` は 1 回だけ、`lean_status` で imports 一致を確認済み。

| 段 | 状態 |
|---|---|
| **段 1**(`j=1`、巡回) | ★**完全に入った**。`exists_natCast_herbrandPhiGroup_of_isCyclic` |
| **段 2**(`j>1` の帰納) | ★**入っていない**(理由は**配管**。下記) |
| **段 3**(`G ≠ G_1`) | ★**入った**。`natCard_quot_dvd_of_herbrandPhi_natCast` |

### ★★段 1 は原文と経路が違う —— しかもそのおかげで退化が自動で通った

★**Corollary 6.7 の p 進展開を経由しなかった。**
`i_{j+1} − i_j` が `p^{j+1}` で割れる(Prop 6.6 (iii) = Y7b)ことから
**直接 `p^m ∣ Σ_{k=1}^n |G_k|`** を出した。
⇒ ★**`n_i ≥ 1` も `1 ≤ j ≤ m−1` も要らず、`m = 0, 1` の退化も自動で通る。**
☆本体は持ち場で「Cor 6.7 を使う」と原文どおり書いたが、**実装者はより短い道を見つけた**
(直近 8 波すべてで同じことが起きている)。

★**Y13 の `dvd_of_sum_eq_natCard_mul` はそのまま差し込めた**(段 3)。
追加は `(G_i).subgroupOf G_1` を `Subgroup.subgroupOfEquivOfLe` で潰す一手だけ。
★**Y12 の `herbrandPhi_natCast` / `herbrandPhiGroup_natCast` もそのまま使えた。**
★**使えなかったのは `herbrandPhiGroup_comp`(Lemma 6.10(ii))だけ** —— 固定環 `C` と
**7 本の適合条件**が要るため、結論の形だけを仮定 `hcomp` として受け取った。

★**思ったより安い**: 見積 500–1000 行に対し **615 行**。効いたのは
(a) 段 1 を「p 進展開の値」ではなく**ブロック和の整除性**に読み替えたこと、
(b) ★**原文が `hence` / `by definition` で畳んだ 3 箇所(ブロック分割・跳びの位置・順分岐の φ)が
いずれも抽象核として切り出せた**こと(**0.08–0.13 秒**)。

### ★★有限アーベル p 群の構造定理は mathlib に在った（ただし未 import）

Y15 が実測(`Unknown constant` で 0.03 秒判別。★#68「無い」ではなく「import していない」):
- `AddCommGroup.equiv_directSum_zmod_of_finite` — `Mathlib/GroupTheory/FiniteAbelian/Basic.lean:135`
- `AddCommGroup.equiv_directSum_zmod_of_finite'` — 同 :151
- 双対: `CommGroup.exists_apply_ne_one_of_hasEnoughRootsOfUnity` — `.../Duality.lean:64`

★**olean はビルド済み**なので **import 1 行で使える**。
★ただし Y15 は「**構造定理があっても段 2 は閉じない**(塔が要る)」と判定して使わなかった。

### ★★★段 2 が止まった理由が、走行中に解消していた

Y15 の判定:
> 段 2 — 入っていない。理由は**配管**(数学ではない): 帰納の各段が
> 「商 `G/H` が固定環 `C` に忠実に作用する」ことを要求し、木の Prop 6.9 / Lemma 6.10(ii) は
> その `C` と **7 本の適合条件**を仮定で受ける形になっている。

★★**その 7 本を、Y16 が Y15 の走行中に供給した**(`fixedRingMulSemiringAction` /
`algebraMap_smul_fixedRing` / `smul_fixedRing_eq_self` ほか)。
さらに Y14 が `fixedRing_isDiscreteValuationRing` / `addVal_map_eq_card_mul_fixedRing` を出している。
⇒ ★**Y15 は Y16 を知らないまま止まった。**

☆★**これは並列化の副作用である**: 2 体を同時に走らせると、
**片方の成果がもう片方の「越えられない」判定を無効化する**ことがある。
★**規約に足すべきこと**: 「越えられない」と報告された持ち場は、
**その波で他に着地したものを突き合わせてから**新ノードにすること。

⇒ **Y15b `HasseArfInduction` として即座に配った**(Y15 と Y16 の両方を import させ、
「Y15 の判定は今は成り立たない可能性が高い。★**ただし本体は実機で確かめていない**。
確かめて、まだ足りなければ具体名で報告すること」と明記)。
★Y17 が `hfix` を**同時に**埋めているので「要るなら仮定として受け取り、
★**Y17 の着地を待たないこと**」も書いた。

`lean-idioms.md` **#118**(`Nat.cast` 経由の `min`/`1` は `min_eq_left` の前に `Nat.cast_one` /
`Order.le_of_lt_add_one` の名前つき引数は `x`,`y` /
`irreducible_iff_uniformizer` は `open IsDiscreteValuationRing` 必須 —— #68 の変種)。

### Y17 着地 —— ★★Λ7 の鎖から最後の仮定が消えた

`Found/PGC/FixedRingMonogenic.lean` **376 行、sorry 0**。
`#print axioms` 5 本すべて `[propext, Classical.choice, Quot.sound]`。
`lake build` 成功(2992 jobs, 6.2 秒)。

主定理 `fixedRing_mem_adjoin_uniformizer` = `hfix` そのもの。
これを Y16 に流し込んだ `card_mul_ramIndex_eq_sum_ramIndex_fixedRing_of_irreducible`
(および `hAne` を自動化した `..._of_injective`)で、**Lemma 6.8 形から `hfix` が消えた**。
⇒ ★★**Λ7 の鎖(Y10 → Y11 Herbrand → Y12 Lemma 6.10 → Y14 → Y16 → Y15 Hasse-Arf)が
原典 §6.1–§6.2 の設定だけに載った。** 今朝は 8 つの仮定を持ち回っていた。

### ★★brief の冒頭警告が、入れた当日に機能した

★**Y17 は `brief.mjs --id lemma-5-11` の冒頭で**
> ★**この項目は既に木にある**: `Found/PGC/LowerRamificationGroup.lean`

**を読み、一般形(`adjoin_uniformizer_eq_top` / `exists_uniformizer_adjoin_eq`)が既にあると
気づいた。** そのうえで
★**「本ノードは同じ原典項目を `K″ = (K′)^H` に当てた実例」**という正しい切り方をし、
`.src` を重ねた(同一 item への複数 `.src` は本木で普通 —— `GenEll p17 Prop 3.4` は
**63 件**と実測して確認したうえで判断している)。
☆**2026-09-07 に本体が踏んだ失敗形([[brief-header-carries-the-warning]])が、
同じ日に道具と持ち場の両方で塞がり、実際に機能した。**

### ★本体の見立てが外れた（今日 9 件目）—— ただし今回も覆った

本体は持ち場に「★`Algebra A C` のインスタンスが無い可能性が高い。
**これが本ノードの重心かもしれない**」と書いた。★**外れた。実質 2 行で立った。**

```lean
@[reducible] def Subring.algebraOfMapsTo (S : Subring B) (h : ∀ a, algebraMap A B a ∈ S) :
    Algebra A ↥S := ((algebraMap A B).codRestrict S h).toAlgebra
haveI : IsScalarTower A ↥S B := IsScalarTower.of_algebraMap_eq fun _ => rfl   -- ★rfl 一発
```
★★**鍵は「大域インスタンスにしない」こと。** `hfix` の statement に `Algebra A C` が
現れないので `letI` で証明の中に閉じ込められる(★#117(i) の害を回避)。
`lean-idioms.md` **#119** に記録。

★**本体の段取りより安かった点が 5 つ**:
(1) `hres` の「向きが違う」問題は**起きなかった**(局所環では非単元の逆像が非単元、**2 行**)
(2) ★**`hadjA`(`𝒪_{K′}=𝒪[π′]`)は §1・§2 で 1 度も使っていない** ——
`K′/K` の完全分岐だけで `K″/K` 版の Lemma 5.11 が出る
(3) ★**`H ⊴ G` は §1・§2 で不要**(固定環は正規性なしに部分環。§3 でだけ復活)
(4) `B` への輸送は **1 行**(`AlgHom.map_adjoin` + `Set.image_singleton` +
`IsScalarTower.coe_toAlgHom'`)
(5) 重かったのは `[Module.Finite A C]` の供給だけで、それも `Module.Finite.of_injective` 1 行

☆逸脱 4 が誠実: ★**`hϖ : Irreducible ϖ` が新たに要る**(Y16 は任意の `ϖ` を取れたが、
`hfix` を供給するには素元でなければならない —— 単元なら `𝒪_{K″} ≠ 𝒪[ϖ]`)。
★**statement をねじ曲げてはいない**(`ϖ` は原典で `K″` の素元)。
逸脱 2・3 で `[IsNoetherian A B]` と `hAne` を足したことも明記
(**仮定を足しただけで弱めていない**)。

☆新ノード候補(急ぎではない): `[IsNoetherian A B]` と `hAne` を Lubin-Tate の具体設定で供給する
——「`𝒪_{K′}` は `𝒪_K` 上有限」。★**数学が足りない側ではなく在庫調べの問題。**

### 新しく配った持ち場: Y18 `Definition 6.12`（上付き番号付け）+ `Corollary 6.13`

★★**pGC 本体に直接効く数少ないノードである。** 構造化係の注記:
> pGC §2(**Definition 2.3、上付き番号付けの高次分岐群**)と**同じ対象を定義**しており、
> ★**Λ7 はこちらを実装本体にして pGC 側をその上に載せる予定。**

★**本ノードの本当の中身は `φ_G` の逆写像の存在**である。構造化係の注記:
> ★**その全単射性は原典では明示されていない**(Lemma 6.10 の形から読み取れ、という扱い)——
> ★**Lean 化では `φ_G` の狭義単調性・連続性・非有界性が独立の節点になる。**

在庫は `strictMono_herbrandPhi`(Y11)だけで、★**連続性と非有界性は在庫に無い見込み**
(★本体は索引で確かめていない)。本体の見当(有限個の `min` の和だから連続 /
`i(1) = ⊤` だから非有界)を**「実機で確かめていない。外れたら訂正して報告すること」**と
明記して配った。

★**本体は今回 `def-6-12` / `cor-6-13` の冒頭を確認してから配った**(どちらも木に無い)。

### Y18 着地 —— 上付き番号付け（Def 6.12）と Cor 6.13 (i)(ii)

`Found/PGC/UpperRamificationGroup.lean` **650 行 / 宣言 41 件、sorry 0**。
`#print axioms` 主要 14 宣言すべて `[propext, Classical.choice, Quot.sound]`。
`lake build` 成功(2990 jobs, 6.7 秒)。`lean_check` 13 回(合計 3 秒弱)＋ `leanfile.mjs` 3 回。
★**`lean_start` も `lean_reset` も呼んでいない**(4 回連続)。

★**本体の見当は 2 つとも当たった**: 連続性は「有限個の `min` の和」で
`Continuous.min` + `continuous_finsetSum`、非有界(上)は `i(1) = ⊤` から `φ_G(n) ≥ −1 + (n+1)/|G|`。
連続性・非有界性は在庫に無く、自分で書いた(4 補題・約 45 行)。

### ★★本体の段取りより 1 段安い道があった（10 度連続の実証）

本体は `StrictMono.orderIsoOfSurjective` / `Continuous.surjOn_Icc` で
`[0,∞)` 上の全単射を作れと書いたが ——
★★**`i(τ) ≥ 1` から `n ≤ 0` では全項が `min{i(τ), n+1} = n+1` になり `φ_G(n) = n`**
(`phiOf_of_nonpos`、**6 行**)。つまり ★**`φ_G` は `(−∞,0]` 上で恒等写像**で下にも非有界。
⇒ ★**`φ_G : ℝ → ℝ` が全域の全単射**になり、`SurjOn`・`Ici` への制限が一切要らず、
`intermediate_value_Icc` **1 回(6 行)**で全射性が出た。逆写像は `Function.invFun` で足りた。

★**(i) は新たな仮定を 1 つも足していない** —— `φ_{G/H}` については単射性(狭義単調性、仮定不要)
しか使わないので `ϖ` が `C` の一意化元であることを足す必要が無かった。
★**商群 `G ⧸ H` は作らずに済んだ(5 波連続)。**

★**(iii) は落とした**(新ノード)。理由: **Theorem 6.11 の一般の可換 `G` に対する形**が要り、
`HasseArf.lean` は巡回群・馴分岐商までしか無い。★**入口の第 1 文だけは入れた**
(`upperRamificationGroup_eq_lowerRamificationGroup`)。
`lean-idioms.md` **#120**(`if` を割った直後の `Continuous fun r => r` を `simp` は閉じない /
`continuous_finset_sum` は deprecated / `a ≤ b → a/c ≤ b/c` は **`gcongr` が一発**)。

### Y15b 着地 —— Hasse-Arf 段 2 が `hind` 以外すべて閉じた

`Found/PGC/HasseArfInduction.lean` **703 行、sorry 0**。
`#print axioms` 全宣言 `[propext, Classical.choice, Quot.sound]`。
`lake build` 成功(3107 jobs, 6.9 秒)。★`lean_start` は呼んでいない(5 回連続)。

★★**本体の見立て「Y15 の『配管が越えられない』は今は成り立たない可能性が高い」は当たった。**
Y15 が挙げた 7 本のうち **`hcomp`・`hHtriv`・`hinj`・`hfixC`・`hres`・`hAC`・
`[MulSemiringAction G C]`・`[IsDiscreteValuationRing C]` が全部消えた**(Y14 + Y16)。

★★**しかし Y15 が気づいていない別の債務があった**(実装者が発見):
段 2 は Y15 の段 1 を**商 `G ⧸ H` が `C` に作用する形**で使うので、`G` の作用(Y16)では足りない。
- `[MulSemiringAction (G ⧸ H) C]` と `hq` → ★**本ノードが払った**(`hq` は `rfl`)
- ★**`FaithfulSMul (G ⧸ H) B^H` → 本ノードが払った**(予定外の収穫。
  Y16 の Lemma 6.8 で `|H|·i_ϖ(σ)=⊤` ⇒ 和の項が `⊤` ⇒ `στ=1`)

### ★構造定理は「インターフェースは巡回商、エンジンは構造定理」に落ち着いた

★段 2 は成分の位数 `p^{m_i}` を **1 度も使わない**。位数は `IsPGroup.to_quotient` で別に出る。
`H := ker ψ` で終わるので分解を経由するより安い。
import は **1 行 `Mathlib.GroupTheory.FiniteAbelian.Basic`**
(★乗法版 `CommGroup.equiv_prod_multiplicative_zmod_of_finite` があるので `Additive` 往復が不要)。
Duality は `HasEnoughRootsOfUnity` のインスタンスが要り高い。

★**Y12 の `herbrandPhiGroup_comp` が今回は使えた**(Y15 では使えなかった)。
差は **`C` を抽象のまま両方の作用付きで持つ**設計にしたこと。

☆逸脱 3 が誠実: **`hind`(`φ_H(n) ∈ ℤ≥0` の再帰呼び出し)は仮定で受ける ——
数学ではなく配管**。
☆予定外に高かった点も正直に: `orderOf σ̄ = p^m` と `σ̄ ∈ (G/H)_1` は**仮定に置かず導いた**
(`h1 : G_1 = ⊤` から)。★**原典「First assume G = G_1」に忠実になった。**

`lean-idioms.md` **#121**(`ℕ∞` に `WithTop.sum_eq_top` を `rw` できない / `mul_top` は無く
`ENat.mul_top` / `Set G` の積に `open Pointwise`)・**#122**(`IsPGroup.of_equiv` はドット記法 /
`Subgroup.zpowers` の所属を開くと β 未簡約の項が出るので `show` を挟む)。

### 新しく配った持ち場 2 つ

**Y19 `RamificationFiltration` の構成**（★**pGC §2 が入力に要求している対象そのもの**）。
`Interface/PGC/LocalFieldData.lean:160` の構造体は
`Gv : (K : PAdicLocalField p) → ℝ → Subgroup K.absGal`(★**絶対 Galois 群の上**)で、
`waiting` に「**Herbrand の定理を要する**」と書いてある。
★★**Y18 の Cor 6.13(i)(`G^m H/H = (G/H)^m`)こそが
「有限 Galois 部分拡大の逆極限として `Γ_K^v` を定義できる」根拠**である。
★**本体が見つけた穴を明示して配った**: Λ7 の鎖は `K′/K` の完全分岐を仮定しているが、
**一般の有限 Galois 拡大は完全分岐ではない**。道 A(惰性群に制限)/ B / C を挙げ、
★**どれが正しいかを測って選ばせる**(本体は確かめていない)。

**Y20 Hasse-Arf 段 2 の残り配管 3 本**。Y15b が名指しした
`Algebra A ↥(fixedRing B H)` / `CharP (ResidueField …) p` / Y17 の `hfix`・`hadjC`。
★**Y17 が `Subring.algebraOfMapsTo` と `adjoin_uniformizer_subring_eq_top` を作り置きしている**
ので、★**そのまま使えるかを測らせる**(本体は型を突き合わせていない)。
★1 を**大域インスタンスにするか `letI` に留めるか**が設計の要
(#117(i) と #119 が衝突しうる)。

### Y20 着地 —— 段 2 の配管 3 本が全部埋まった（337 行、思ったより安い）

`Found/PGC/FixedRingBaseAlgebra.lean` **337 行、sorry 0**。
`lake build` 成功(3109 jobs, 9.4 秒、警告 0)。`leanfile.mjs` **4 往復・実質 1 エラーのみ**。
★抽象核 `smulCommClass_of_smul_algebraMap` は **`does not depend on any axioms`**。

| Y15b が名指しした債務 | 供給した宣言 |
|---|---|
| `Algebra A ↥(fixedRing B H)` | `fixedRingAlgebra` |
| `SMulCommClass (G ⧸ H) A C` | `smulCommClass_quotient_fixedRing` |
| `hresC` | `exists_sub_mem_maximalIdeal_fixedRing` |
| `CharP (ResidueField C) p` | `charP_residueField_fixedRing` |
| `hfix` / `hadjC` | ★**Y17 の 2 本がそのまま嵌った** |

★★**段 2 の仮定は「原典 §6.1–§6.2 の底の設定」と `hind` だけになった。**

### ★★本体の見立てが外れた（今日 11 件目）—— 今回は「要ると書いた仮定が要らなかった」

本体は「2（`CharP`）は `K′/K` が完全分岐なら剰余体が同じ」と書いたが ——
★★**`hresA`（完全分岐）は 2 には要らなかった。** `IsLocalHom` すら要らず、
Y17 の `mem_maximalIdeal_of_map_mem`（単元の像は単元、という**易しい向き**）と
`CharP.charP_iff_prime_eq_zero` の **2 本で 3 行**。
★**不分岐でも剰余体が伸びるだけで標数は `p` のまま。**
⇒ `charP_residueField_fixedRing` は `hresA` を受け取っていない。

### ★`Algebra A C` は `letI` に留めた —— 測って決めた

本体は「大域インスタンスにするか `letI` に留めるかが設計の要」と書いた。実装者の判定:
> 系の**結論**に `Algebra A C` が現れず、現れるのは `hadjC`/`hresC` 側だけで、
> それは本ファイルが自分で供給する ⇒ **系の statement から `Algebra A C` が完全に消える**。
> ★**逆に大域にすると `A := ℤ` で `Int.instAlgebra` とダイヤモンドになる。**
> ★**`#117(i)` とは衝突しなかった** —— `fixedRing B H : Subring B` は型としては
> `↥(Subring …)` なので `Subring.toAlgebra`・`IsDomain ↥S` は普通に降りてくる
> （`fixedRing` を**定義まで開く場面が無い**）。

★**Y17 の `Subring.algebraOfMapsTo` と `adjoin_fixedRing_uniformizer_eq_top` は
新しく書かずにそのまま使えた**（`letI` の入れ方が Y17 と同一項なので
`exact` の既定透明度で一致）。

★**`lean_start` を呼ばずに済ませたのは 6 回連続。** Y20 の判断が良い前例:
自分の在庫が共有 env の**部分集合でない**と判明したので、
**mathlib のみの `#check`・抽象核だけ共有 REPL で測り、ABC3 在庫を引く部分は
`leanfile.mjs` に回した**。Y19 の環境は無傷。

`lean-idioms.md` **#123**: ★**`haveI := f (A := A) (fun σ c => …)` のように
引数がラムダだけだと `B`・`G`・`H` が決まらず `typeclass instance problem is stuck` が出る**
——`haveI :=` には**期待型が無い**のが原因。暗黙引数を名前で全部渡せば直る
（同じ項を `exact` の引数位置に置くと期待型があるので何も渡さずに通る）。
併せて「局所環の代数写像に沿って剰余体の標数は降りる（`IsLocalHom` 不要・分岐不要）」を在庫として記録。

### 新しく配った持ち場: Y21 `hind` の除去 —— ★Hasse-Arf が閉じる最後の 1 本

Y15b と Y20 の実装者が**揃って**言っている:
> `hind` を消すには「**部分群 `H ⊆ G` すべてについて `φ_H` の主張を強帰納法で回す**」
> 新ノードが要る（各段で塔 `K′/K^H/K` を立て直すため）。
> ★**Y20 の 3 本が入ったので、その新ノードは底の設定だけから塔を組めるようになった。**

★本体の見当（★**実機で確かめていない**と明記して配った）:
- 帰納の測度は **`Nat.card H` の強帰納法**（`H ⊊ G` ⟹ `|H| < |G|`）
  ——原文は「`j`（巡回成分の個数）」だが、★**Y15b は成分の位数を 1 度も使っていない**
- `H = H_1` は `H ≤ G = G_1` から出るはず

★**閉じれば Theorem 6.11（Hasse-Arf）が 3 段そろって完成する。**

---

## ★★★Λ6 §4 を `L = K̂` で立てる道が測定で決着した（読み取り専用の測定、2026-09-07）

本体は「Λ6 §4 の実体は『新しい数学』ではなく『`L = K` から `L = K̂` への一般化』かもしれないが、
★**規模を測っていない**」として在庫調査係に測らせた。**答えが出た。**

### 測定 1: Lubin-Tate 在庫 **57 ファイル / 14,525 行 / 471 宣言**の担い手

| 分類 | 件数 | ファイル | 行数 |
|---|---|---|---|
| `PAdicLocalField p` 特化 | **194** | 26 | 7,301 |
| **抽象環**(`{A}[CommRing][IsLocalRing][IsDomain]…`) | **269** | 31 | 7,224 |
| 担い手なしの汎用補題 | 8 | — | — |

★**主要 4 宣言のうち 3 つは抽象環**(`iteratedLubinTateDistinguished` / `LubinTateEndo` /
`formalGroupLaw`)。`iteratedLubinTateTorsionPoints` が `K` 特化なのは
**根を取る先に代数閉体 `K.closure` が要るから**だけ。

★測定の工夫として記録に値する: **「statement + その宣言より前にある `variable` 行」を
突き合わせないと 121 件が未分類になる**(`variable` 束縛は索引の statement 欄に現れない)。

### ★★★測定 2: `unramifiedCompletion K` は `PAdicLocalField p` に **ならない**（2 つの独立な理由）

**理由 1（定義から）** `Skeleton/PGC/Setup.lean:40` の `PAdicLocalField` は
`[FiniteDimensional ℚ_[p] carrier]` を要求する。`K̂^ur` は `K^ur` を含むので **ℚ_p 上無限次元**。

**理由 2（型で矛盾が出る）** ——★**これが決定的**:
- Lubin-Tate 在庫の抽象環 269 件のうち **100 件**が
  `[Fintype (IsLocalRing.ResidueField A)]` + `hq : card (ResidueField A) = pp^ff` を要求する
- しかし `isAlgClosed_residueField_unramifiedCompletionInt`(`UnramifiedResidueField.lean:935`)
  により **`𝒪_{K̂^ur}` の剰余体は代数閉**
- mathlib の `instance [Field K] [IsAlgClosed K] : Infinite K` と衝突する

★**`lean_check` で実測(0.01 秒)**:
```lean
example (κ : Type) [Field κ] [IsAlgClosed κ] [Fintype κ] : False := by
  have := (inferInstance : Infinite κ); exact not_finite κ
```
⇒ ★★**「既存の Lubin-Tate 機械を `K̂` に代入するだけ」は論理的に成立しない。**

☆未測定として明記された点: `𝒪[unramifiedCompletion K]`(ノルムの付値環)と
`unramifiedCompletionInt K` を同一視する定理は木に **0 件**。理由 2 を Lean で完結させるなら 1 本要る。

### ★★測定 3: 有限剰余体が**本質的に**効いているのは 2 箇所しかない

`hq` の**本質的**使用(引数の受け渡しではない)は corpus 全体で **14 箇所**:
`1 < pp^ff` に使うもの **12 件** / `0 < pp^ff` に使うもの **2 件**。
★**この 14 箇所は `hq1 : 1 < pp^ff` に置き換えれば有限性が不要になる。**

★**真に有限体を使うのは 1 本だけ**: `mvPowerSeries_pow_card_eq_expand`
(`Found/PGC/FrobeniusExpand.lean:47`、**59 行**)と、その呼び出し側 2 箇所。
これは原典 Lemma 3.4 の合同を **φ = id と特殊化**したもの。

★★★**その φ 捻り版は既に木にある** —— Λ5b/Λ6 が同じ壁に当たって作っていた:
`pow_subst_eq_subst_map_iterateFrobenius`(`DworkThetaStep2.lean:448`)。
★**有限性を仮定しない**(任意の `[CommRing R] [ExpChar R pp]`)。
★`DworkThetaStep2.lean:81` と `:446` のコメントが
**「`mvPowerSeries_pow_card_eq_expand` は有限体専用なのでここでは使えない(`𝓀_{K̂^ur}` は無限)」**と
明示的に書いている。☆**過去の自分が既に同じ壁に当たり、回避策を残していた。**

### ★★★決定（D28）: **道 A（基底変換）を採る**

| | 道 A（基底変換） | 道 B（原典の一般性 `f ∈ 𝒪_L[X]`） |
|---|---|---|
| 中身 | `f ∈ 𝒪_K[[X]]` のまま `L^m_f := K̂^ur(µ_{f,m})` を `ℂ_K` の中で取る | φ 捻りを 269 宣言の**型**に通す |
| 既存への影響 | ★★**Lubin-Tate 機械 471 宣言を 1 件も触らない** | 広い署名変更（件数は**未測定**） |
| 木の現状 | ★**既にこの道を採っている**(`lubinTateCompletionField`、`LubinTateTowerFIndependent.lean:364`) | — |
| 数学的な壁 | 無し | 1 本(`FrobeniusExpand`)。★**代替は既に在る** |

★★**道 A を採る。理由**: (1) 既存 14,525 行に一切触らない (2) ★**木は #9(Cor 4.9 前半)で
既にこの道を採っており、`lubinTateCompletionField` がその実体である** (3) 道 B は
数学的な壁こそ無いが**署名変更の件数が測れていない**。
☆**これは技術的な経路選択であり人を待つ判断ではない**ので、決定として記録する。

### ★道 A で要るのは 3 本だけ（測定で確定）

1. **`[K̂^m_f : K̂^ur] = q^n − q^{n−1}`** —— 木に **0 件**。
   材料は `totallyRamifiedAdjoin_inf_unramifiedClosure`(`TotallyRamified.lean:687`)と
   `inf_eq_bot_of_isUnramified_of_isTotallyRamified`(`RamifiedUnramifiedDisjoint.lean:46`)にあるが、
   ★**`K^ur → K̂^ur` の完備化を渡る段が無い**
2. **`Gal(K̂^m_f/K̂^ur) ≃* (𝒪_K/π^n)^×`**（＝ Prop 4.4(iii) の `L = K̂` 版）—— 木に **0 件**
3. **`N(−α) = π`** —— 木に **0 件**。
   ★`Algebra.norm` を型に含む宣言は木全体で **4 件**(Arakelov 3・SixExp 1)、`Found/PGC` には **0 件**

☆★**在庫調査係は「私は見積もっていない」と明記した**(`L=K` 版の対応物が
68 + 153 + 425 行であることは書いたが、それは類似物の行数であって見積もりではない、と)。
★**この誠実さは正しい。** 本体も見積もりを捏造しない。

### ★その他の測定結果

- **mathlib に Lubin-Tate は無い**(`grep -ci lubin` → **0**、
  `absent-recheck.mjs --try 'lubin.?tate'` → **0 件**)。
  最も近いのは `FormalGroup`(`Mathlib/RingTheory/FormalGroup/Basic.lean:53`)の **12 宣言のみ**で、
  形式 `O`-加群も Lubin-Tate 級数も無い。
- **`K̂` 側の在庫は 149 件 / 11 ファイル**。`𝒪_{K̂^ur}` の環論的構造は既に揃っている
  (`isDiscreteValuationRing_` / `isAdicComplete_` / `maximalIdeal_…_eq_span` /
  `expChar_residueField_` / `isAlgClosed_residueField_`)。
- **Prop 4.4(i) の `𝒪/𝔭^m ≃ₗ µ_{f,m}` という 1 本の同型は木に無い**(部品は揃っている)。
- ★**Prop 4.7(ii) の原文は「Let L = K̂」で確定**(`.txt` 行 421 が `(ii) Let L = bK.`)。

### Y19 着地（部分）—— `RamificationFiltration` の抽象核が完全に埋まった

`Found/PGC/AbsGalRamificationFiltration.lean` **577 行、sorry 0**（見積 600–1200 に対し ★安い）。
`#print axioms` 全 20 宣言が `[propext, Classical.choice, Quot.sound]`
（`comap_compat_of_coe_mul_coe_eq` と `directed_openNormalBase` は **`Classical.choice` すら不要**）。
`lake build` 成功(3021 jobs, 6.4 秒)。★`lean_start` は呼んでいない（7 回連続）——
共有 REPL で `lean_check` 9 回（合計 **1.5 秒弱**）、`leanfile.mjs` 7 往復。

★★**`ramificationFiltrationOfStages : (∀ K, StageFiltration K.absGal) → RamificationFiltration p`
が出て、4 フィールド（`Gv`/`isClosed`/`isNormal`/`antitone`）全部が埋まった。**
⇒ ★**あとは各有限段のデータを 1 点作れば、pGC §2 の入力の本物が出る。**

★**抽象核は位相群の言葉だけ**（分岐・付値・Galois が 1 語も出ない）で **0.03–0.36 秒**。
最大の山と目された**コンパクト性論法**（`coe_limit_mul_coe`：極限は各段へ全射）が
★**1 往復 + 1 修正**で通った。一意性 `eq_limit_of_isClosed`（各段へ全射する閉部分群は極限に限る）は
**コンパクト性すら不要**。

### ★★★本体が挙げた「道 C」を、実装者が**反例で否定した**（記録に値する）

本体は完全分岐の穴の埋め方として道 A / B / C を挙げた。Y19 の判定:

> **道 C は使えないと確定**: ★**完全分岐な有限 Galois 部分拡大の族は有向でない。**
> 反例 `ℚ_3(√3) · ℚ_3(√−3) ∋ √−1`、`ℚ_3(√−1)` は**不分岐**
> （`−1` は mod 3 で平方非剰余）。**逆極限の底にできない。**

☆★**本体は 3 つの道を挙げただけで、どれが偽かを知らなかった。**
**具体的な反例で潰されたのは今日初めて**であり、
「測って選ばせる」形式の価値がまた 1 つ実証された。

⇒ Y19 は **道 A/B（惰性部分に寄せる）**を選び、`L₀ := L ∩ K^ur` として
`L/L₀` に Yoshida §6 の設定を当てる方針にした。

### ★Y18 の結論の「形」がそのまま効いた（設計の連鎖）

> `upperRamification_coe_mul_coe_eq` の**結論の形がそのまま `compat` の形**だった
> ——★**Y18 が商群を作らず集合積で書いたのが効いた。**

☆Y9・Y10・Y11・Y12・Y18 と **5 波続けて商群を作らなかった**選択が、
**6 波目（Y19）の逆極限の定義をそのまま可能にした。**
★ただし**実際には代入していない**——仮定 10 個が各有限段で完全分岐を要求するため、
橋（`comap_compat_of_coe_mul_coe_eq`）だけ先に置いた。

★**退化 witness を作って明示した**: `trivialRamificationFiltration`（`v>0` で `⊥`）。
★**「原典の `Γ_K^v` ではない」「G2 に使うな」と docstring に書いてある**
——☆**構造体を満たす項が存在することを示しつつ、それが本物でないと自分で釘を刺した。**

`lean-idioms.md` **#124**: ★**`AlgEquiv.restrictNormalHom` と `restrictNormalHom_surjective` は
`E` と `K₁` の役割が逆**（素直に名前付き引数を書くと `failed to synthesize Algebra E ↥L`）/
`∃ P, (IsOpen (P : Set Γ) ∧ P.Normal)` は型注釈が無いと `P : Γ → Prop` に潰れ
**`Function.Normal` という謎エラー**になる。

### 新しく配った持ち場: Y19b+c 惰性への還元

Y19 が名指しした残りの数学 2 つ:
> (a) 下付き分岐群が**不分岐底変換で不変** ——★実装者いわく「**付値が同じという理由だけ**」
> (b) `Gal(L′/L′₀) ↠ Gal(L/L₀)` の全射性と、その商が Cor 6.13(i) の `G/H` に一致すること

★Y19 が **`surjective_restrictNormalHom`** と **`ker_restrictNormalHom_eq_fixingSubgroup`**
（体論の部分）を既に作ってあるので、**`L₀` 側に移すのが仕事**のはず。
★本体の見当「`G_i ⊆ G_0 = Gal(L/L₀)` だから一致する」は
**実機で確かめていない**と明記して配った。
★§3（段データの組み立て = Y19d）は余力扱い。
★★**これが済めば `Skeleton/PGC/Section2` の入力の本物が出る。**

---

## ★★★★★Yoshida 2008 の証明に**論理の穴**を見つけた（erratum の 3 件目。今度は誤植ではない）

Y21（`HasseArfStrongInduction.lean`、491 行、sorry 0）が段 2 を閉じる過程で発見した。

### 原文（Theorem 6.11、`.txt` 1102–1124 行）

> For j > 1, if `G_n ≠ G_{n+1}` we can find H with `G/H ≅ ℤ/p^{m_i}ℤ`, and `G_nH/H ≠ G_{n+1}H/H`.
> **We have `φ_H(n) ∈ ℤ≥0` by inductive hypothesis**, …

★★**この「by inductive hypothesis」が引けないことがある。**
帰納法の仮定は **`H_n ≠ H_{n+1}`** を要求するが、`H_i = H ∩ G_i` なので**成り立たない場合がある**。

### ★実装者が挙げた反例

`G ≅ (ℤ/p)²`、`G = G_1`、跳びが 2 つ `n₁ < n₂`（`G_i = C`（位数 p）for `n₁ < i ≤ n₂`）。
`n = n₂` のとき原文の `H` は `|H| = p` かつ `C ⊄ H` ⇒ `H ≠ C` ⇒ `H ∩ C = 1`
⇒ ★**`H_{n₂} = H_{n₂+1} = 1`**（＝帰納法の仮定が引けない）。

★★**しかも `φ_H(n₂) = n₁ + (n₂−n₁)/p` で、その整数性は結論そのものと同値である。**
⇒ ★★**原文はこの場合、結論と同値のものを「帰納法の仮定」として引いている。**

### ★★★そして実装者は帰納法ごと消して穴を塞いだ

`Nat.strong_induction_on` は **1 つも書いていない**。
`exists_natCast_herbrandPhi_of_jump` が **Proposition 6.9 だけから** `φ_H(n) ∈ ℤ` を出す:

> `t := i_ϖ(x) − 1` は**整数**で `φ_H(n) ≤ t < φ_H(n+1)`。ここで
> `|H|·φ_H(n) = Σ_{i=1}^n |H_i|`、`|H|(φ_H(n+1) − φ_H(n)) = |H_{n+1}|` であり、
> `H_{n+1} ≤ H_i` (i ≤ n) と Lagrange から **`|H_{n+1}|` は `Σ|H_i|` と `|H|` の両方を割る**。
> ★**長さ `|H_{n+1}|` の半開区間に入る `|H_{n+1}|` の倍数は左端だけ** ⟹ `φ_H(n) = t`。

★**抽象核は純算術 2 本 + 純群論 1 本**（`eq_of_le_of_lt_add_of_dvd` /
`div_eq_natCast_of_dvd_of_le_of_lt` / `natCard_dvd_sum_natCard`）で、共有 REPL **0.12 秒**。
`eq_of_le_of_lt_add_of_dvd` は `#print axioms` が **`[propext, Quot.sound]` のみ**。

★**見積 400–800 行 → 491 行だが内訳が全く違う**（強帰納法 **0 行**、代わりに整除性 3 本）。
往復は合計 **6 回**（うち一発が 2 回）。★**本体の段取り（`Nat.card` の強帰納法）は
「測度自体が不要になった」という形で外れた**（今日 12 件目）。

☆★**この論文で見つけた erratum は 3 件目**だが、前 2 件（Lemma 6.5 の 1 ずれ /
Prop 6.6 の "Let j = 1"）は**誤植**であり、★**今回は論理の穴である**。
★**Lean 化がなければ気づけなかった種類の発見**である。
★**ユーザーに報告すべき事項として記録する。**

### Y21 のその他

★**`brief.mjs` の 1b が発見の入口だった**（実装者の証言）:
> 冒頭で「既に木にある」が分かり、**1b の `.txt` 抽出（`̸=` が生きている形）**で
> `H_n ≠ H_{n+1}` が帰納の前提であることが読め、★**それが穴の発見につながった。**

☆**`pdftotext` は `≠` を `=` に潰す**（[[pdftotext-drops-negation]]）。
★**PyMuPDF 製の `.txt` は `̸=` を残す**ので 1b では読めた ——
**抽出器の違いが数学的発見の分かれ目になった実例。**

`lean-idioms.md` **#125**（★`SetLike` の台集合等式を `∈` ゴールに `rw` できない /
`Nat.dvd_sub'` の消滅 / `MulSemiringAction G ↥(fixedRing B H)` に `[H.Normal]` が要る）。

### 新しく配った持ち場: Y22「部分群 `H` を群として見る橋」

★**Hasse-Arf が 3 段そろう最後の 1 本**。段 3（Y15）が持つ `hk : φ_{G_1}(n) ∈ ℤ≥0` を
段 1+2 から供給するために要る。Y21 が名指しした 4 つ:
`MulSemiringAction ↥H B`（**mathlib にインスタンス無し**と Y21 が報告）/
`herbrandPhi α H n = herbrandPhiGroup ↥H α n` /
`lowerRamificationGroup B ↥H i = (… G i).subgroupOf H` /
`FaithfulSMul ↥H B`・`SMulCommClass ↥H A B` の移送。

★**手がかりとして `herbrandPhiGroup_eq_herbrandPhi_top`（`⊤` の場合の橋）が既にある**ことを
持ち場に書いた（★**本体の見当。実機で確かめていない**と明記）。
★Y12 の実測「`Fintype ↥(⊤ : Subgroup G)` は `[Fintype G]` から推論されない」も警告した。

### Y19b+c 着地 —— 惰性への還元（632 行 / 34 宣言、sorry 0）

`Found/PGC/InertiaReduction.lean`。`lake build` 成功(3083 jobs, 8.9 秒)。
`#print axioms` 全 20 宣言が `[propext, Classical.choice, Quot.sound]`。
★`lean_start` は **1 回も呼んでいない**(9 回連続)——共有 REPL の `lean_check` **18 回** +
`leanfile.mjs` 8 往復。

★**§1 の抽象核は原典より仮定が弱い形で立った**: `𝔪.inertia G ≤ H` なる**任意の** `H` で
`Subgroup.map H.subtype (G_n(H)) = G_n(G)`。**分岐・付値・Galois の語彙ゼロ、0.08 秒。**
§2 の抽象核(純群論)は **0.02 / 0.45 秒・一発**。

### ★★本体の段取りが「必要より強い要求」だった（今日 13 件目の訂正）

本体の見当の当否（実装者の報告）:
- 「`G_i ⊆ G_0` だから `G_0 ≤ H` で一致」→ ★**当たり**(`subgroupOf` が `rfl` なので 3 行)
- 「`G_0 = 惰性群`」→ ★**当たり**(`lowerRamificationGroupAdjoin_zero_eq_ker_residueGalHom`)
- 「`惰性群 = Gal(L/L₀)`」→ ★**半分だけ当たり**。証明できたのは **`Gal(L/L₀) ≤ G_0` の片側だけ**

★★**しかし消費側が要るのは半分のほうだった**:
`lowerRamificationGroup_inertiaGal_zero_eq_top : G_0(L/L₀) = ⊤`、
すなわち **`L/L₀` は完全分岐**が `Gal(L/L₀) ≤ G_0` **だけ**から出る。
⇒ ★★**本体の段取り (a)「`G_i` が不分岐底変換で不変」は必要より強い要求だった。**
☆**「外れた」ではなく「過剰だった」という訂正であり、実装者が消費側から逆算して気づいた。**

### ★Y19 の在庫を使わずに済ませた判断

`surjective_restrictNormalHom` / `ker_restrictNormalHom_eq_fixingSubgroup`(Y19 が作ったもの)は
**使わなかった**。両方 mathlib に同等物があり
(`AlgEquiv.restrictNormalHom_surjective` / `IntermediateField.restrictNormalHom_ker`)、
`L₀` 側への移送は `Subgroup.comap_map_eq` + `map_sup_of_le_ker` の純群論で済んだ。
★★**中間体の塔(`↥L` 上の `↥L'`)を一切作っていない**（`lean-idioms.md` **#59** 回避）。
☆**同じ木の中の新しい在庫より mathlib の方が軽かった、という珍しい例。**

### ★★退化 witness が 2 つ揃った（どちらも「本物ではない」と明記されている）

- Y19 の `trivialRamificationFiltration`（`v>0` で `⊥`）＝ ★**下からの退化**
- Y19b+c の `inertiaStageFiltration`（`v ≥ 0` で `I_K`）＝ ★**上からの近似**

★★**両方の実装者が自分で docstring に「原典の `Γ_K^v` ではない」「本物として使うな」と
書いている。** ☆**構造体を満たす項を作ってしまえる場所で、自分から釘を刺す作法が
2 波続けて守られた。**

### ★実装者が正直に申告した残りの穴 2 つ

1. `G_0 = Gal(L/L₀)` の**逆包含**（＝等式）は未証明。
   ★消費側が要るのは片側だけなので**急ぎではない**、と実装者自身が判定している。
   埋めるには還元射 `Gal(L/K) → Gal(k_L/k)` の**全射性**が要り、在庫の
   `residueGalHom_bijective` は `IsUnramifiedAdjoin` を仮定するので一般の `L` に当たらない。
2. 合成体の Galois 群の同型 2 つ（`Gal(L'/L'₀) ⊓ Gal(L'/L) = Gal(L'/L·L'₀)`）は未証明。

`lean-idioms.md` **#126**（★`IntermediateField.fixingSubgroup` の `Normal` は
**mathlib のインスタンスではない** —— 商群が**宣言の型**に出ると `haveI` では間に合わないので
**前もって `instance` を置く**）・**#127**（★`iff` の前に明示引数があると `.mp` が
「Unknown constant」になる。`rw` は通るのに term モードで落ちる）。

### 新しく配った持ち場: Y19d 段データの組み立て —— ★pGC §2 の入力の**本物**

★3 波（Y18 → Y19 → Y19b+c）でここまで来た。★**残るのは「各有限段のデータを 1 点作る」だけ。**
★★**持ち場に「近似で埋めて『できた』と言わないこと」を明記した** ——
退化 witness が 2 つ既にあるので、**本物と取り違える危険が実在する**ため。
★Y19 が「仮定 10 個が各有限段で完全分岐を要求するので代入できていない」と書いた障害が
**Y19b+c の `L/L₀` 完全分岐性で消えたか**を確かめさせる（★本体は確かめていない）。

### Y22 着地 —— ★★Hasse-Arf の 3 段が合流した（419 行、うち証明は 30 行弱）

`Found/PGC/SubgroupActionBridge.lean` **419 行、sorry 0**。
`lake build` 成功(3111 jobs, 6.2 秒)。`leanfile.mjs` **4 往復・本体は一発**。
★`lean_start` / `lean_reset` は **1 度も呼んでいない**(10 回連続)。
`#print axioms` 全宣言 `[propext, Classical.choice, Quot.sound]`(`comm_subtype` は `[propext]` のみ)。

### ★★★Y21 の「mathlib に無い」は誤りだった —— 4 つとも在った

Y21 は 4 つを名指しして「`MulSemiringAction ↥H B` は **mathlib にインスタンス無し**」と報告した。
★**Y22 が測ったら 4 つとも在った**:

| 要るもの | 実際 |
|---|---|
| `MulSemiringAction ↥H B` | ★`Subgroup.mulSemiringAction`(`Mathlib/Algebra/Ring/Action/Subobjects.lean:40`)。`exact?` が即答 |
| `SMulCommClass ↥H A B` | ★`Subgroup.smulCommClass_left`(`Algebra/Group/Subgroup/Actions.lean:43`) |
| `FaithfulSMul ↥H B` | ★無名 instance(同 `:59`) |
| `lowerRamificationGroup` の一致 | ★`AddSubgroup.subgroupOf_inertia`(`Algebra/Group/Subgroup/Basic.lean:1077`)。`Ideal.inertia` が reducible なので**定義展開すら不要** |

★**引き当ては名前空間 1 回 grep**(`grep -n "FaithfulSMul" .cache/mathlib-index.txt | grep -i subgroup`)。
★**新しい `instance` は 1 つも書いていない**(2 本立ちを避けるため)。
☆★**「無い」の報告は実装者からでも疑うべき**、という実例。
本体は Y22 の持ち場に「★**あなた自身で確かめること**」と書いており、それが効いた。
⇒ ★**Y23 の持ち場に「mathlib を先に引いたか」を返す項目として入れた。**

### ★★本体の見当が外れた（今日 14 件目）—— しかも向きが逆だった

本体は「`herbrandPhiGroup_eq_herbrandPhi_top`(`⊤` の場合の橋)が既にあるから、
**一般の `H` でも同じ道が通るはず**」と書いた。★**外れた。一般の方が易しかった。**

```lean
  rw [herbrandPhi_eq_phiOf, herbrandPhiGroup]
  rfl            -- ★これで閉じる
```
★★**`⊤` の困難は「一般化が難しい」のではなく `↥⊤ ≠ G` という型の食い違いだった。**
`⊤` 版は `Fintype.sum_equiv Subgroup.topEquiv` + `Nat.card_congr` を要していた。
`lean-idioms.md` **#129** に「★`⊤` 版が既にある補題は**一般の方が易しいことがある。
まず `rfl` を叩く**」として記録された。

★Y12 が警告した `Fintype ↥(⊤ : Subgroup G)` の穴には**落ちなかった**
(`[Fintype ↥(lowerRamificationGroup B G 1)]` を段 3 と同じくインスタンス引数として受け取ったため)。

### ★★思ったより桁違いに安い —— 4 項目のうち実質の証明は 1 つだけ

見積 300–600 行に対し **419 行だが、うち証明は 30 行弱**で残りは docstring。
★**実質的な証明を要したのは `hne` の移送(`subgroupOf_injOn` + antitone、**5 行**)だけ。**
その分を余力扱いだった §5・§6・§7 に回し、★**3 段の合流
`exists_natCast_herbrandPhiGroup_of_abelian` まで作った。**

☆逸脱 1 が誠実: ★**`hne` の移送は `1 ≤ n` でしか成り立たない**(数学であって配管ではない)。
`n = 0` では `G_0` と `G_1` の `H = G_1` への制限がどちらも `⊤` に潰れる。
★**原文も別扱いにしている**("If n = 0 then φH(0) = 0.")ので §4 は場合分けを持つ。
☆★**その 1 文を `brief.mjs` の 1b が教えた** ——実装者いわく
「★**これが無ければ `n = 0` で `hne` の移送に失敗して数往復していた**」。

### Hasse-Arf の現状 —— ★段 1+2 は完全に閉じた

- 段 1+2 側は**完全に閉じた**(`G = G_1` 分岐では `C` も `hC` も一切使わない)
- 段 3 側は残る入力が **`hcomp` / `honeC` の 2 つだけ**で、
  ★**`G ≠ G_1` のときにだけ要求される形**に絞れた
- 3 つ目の `htopC` は Y22 が供給済み(`ramIndex_fixedRing_eq_top_of_mem`、
  ★**`ramIndex_eq_top_iff` + `mem_fixedRing` を 1 回ずつ、分岐の議論ゼロ**)

`lean-idioms.md` **#128**(部分群への作用の制限は mathlib に全部ある)・
**#129**(上記)・**#130**(`rw [Nat.cast_zero, Nat.cast_zero]` は 2 つ目で必ず落ちる)。

### 新しく配った持ち場: Y23「順分岐商 `K′^{G_1}/K` の塔」

★**Hasse-Arf が仮定ゼロで閉じる最後の 2 本**(`hcomp` / `honeC`)。
★Y22 が「`htopC` は供給済み、`C := fixedRing B G_1` と `[H.Normal]` まで**道が付いている**」と
書いているので、そのまま配った。
★`honeC`(`σ ∉ G_1 ⟹ i_ϖ(σ) = 1`)は **`htopC` と同じ道が通るか**を測らせる。
★★**「mathlib を先に引いたか」を返す項目に入れた**(Y21 の誤報告の直後なので)。

---

## ★★★★★Theorem 6.11（Hasse-Arf）が**仮定ゼロで閉じた**（2026-09-07、4 波）

`Found/PGC/TameQuotientTower.lean`（Y23）**445 行・10 定理、sorry 0**。
`#print axioms` 6 本すべて `[propext, Classical.choice, Quot.sound]`。
`lake build` 成功(3112 jobs, 6.3 秒)。★`lean_start` は呼んでいない（**11 回連続**）。

```lean
exists_natCast_herbrandPhiGroup_of_abelian_of_setup
    … (hπ' : Irreducible π') (hresA) (hadj) (hAinj)
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    ∃ k : ℕ, herbrandPhiGroup G π' (n : ℝ) = (k : ℝ)
```
★★**残る仮定は原典 §6.1–§6.2 の底の設定と `G` 可換だけ。`hC` は消えた。**
非負性は `herbrandPhiGroup_nonneg_of_abelian_of_setup`。

### ★4 波の内訳（1 日で組み上がった）

| 波 | 行 | 出したもの |
|---|---|---|
| **Y15** | 615 | 段 1（巡回）・段 3（`G ≠ G_1`）。★段 1 は **Cor 6.7 を経由しない**短い道 |
| **Y21** | 491 | ★段 2 を **`hind` ごと**（**帰納法不要**）。★★**原典の論理の穴を発見** |
| **Y22** | 419 | ★段 3 の `hk` を段 1+2 から供給。3 段の合流。★**証明は 30 行弱** |
| **Y23** | 445 | ★`hcomp` / `honeC`。★**仮定ゼロで閉じた** |

★★**Hasse-Arf は Milne も Serre 1961 も Sharifi も証明を外注していた唯一の定理**である
（2026-09-06 の文献調査で確定していた）。

### ★本体の見立ての当否（Y23）

- 「`hcomp` は新ノード扱いになるかも」→ ★**遥かに安かった**。
  **10 個の適合条件は Y14・Y16・Y17 の在庫でそのまま供給でき**、
  `herbrandPhiGroup_comp` への**1 回の代入（項 1 本、証明 5 行）**。
  `hϖ : ι ϖ = π″` は `rfl`、`hcomp`（包含の同変性）も Y16 が `rfl` で出していた。
- 「`honeC` は付値の定義に戻るだけ / 分岐の議論は要らない」→ ★**外れた**（今日 15 件目）。
  `htopC`（Y22）は `ramIndex_eq_top_iff` + `mem_fixedRing` で済んだが、
  ★**`honeC` は `(G/G_1)_1` が p 群であること（`isPGroup_lowerRamificationGroup_one`）を
  経由しないと出ない。** ただし経由さえすれば残りは 3 行。

★**`brief.mjs` の 1b が決定打だった**（実装者の証言）:
> `As φG/H(n) = n/e0 for n ∈R≥0 by definition` の `by definition` が
> 「順分岐商の第 1 分岐群が自明」に等しいこと、そして直後の
> `As e0 and |H| are coprime` が **§1 の抽象核そのもの**であることを教えてくれた。
> ★**これが無ければ `φ_{G/H}(n)=n/e_0` を実数の等式として証明しようとして数往復していた。**

★**mathlib を先に引き、新しい `instance` / `def` を 1 つも書いていない**
（Y21 の誤報告 → Y22 の訂正の直後に配った持ち場で、指示が守られた）。

`lean-idioms.md` **#131**（`Eq.ge` の `.le` は落ちる）・**#132**（`ENat.one_le_iff_ne_zero` 非推奨）・
★★**#133**（**`letI`/`haveI` は宣言の境界を越えない** —— 抽象核と具体層に割ると
`haveI` で入れたインスタンスが片側にしか残らない。
★**Y23 の唯一の失敗がこれで、数学ではなく分割の副作用だった**）。

☆**再利用のために独立させた 2 本**（実装者の配慮）:
`lowerRamificationGroup_quotient_fixedRing_one_eq_bot`（順分岐商の第 1 分岐群が自明）と
★`herbrandPhiGroup_comp_fixedRing`（**`H` が `G_1` である必要が無く任意の有限正規部分群で成り立つ形**）。

### 新しく配った持ち場: Y24 `Corollary 6.13 (iii)`

★**Y18 が「一般の可換 `G` に対する Hasse-Arf が要る」として落とした項目**。
★**その前提がいま満たされた。**
★本体の見当（θ_0 で `q−1` / θ_n で `q` / 跳びは高々 `m−1` 回 / 望遠鏡積）を
**「実機で確かめていない。外れたら訂正して報告すること」**と明記して配った。
★`m = 0` の退化（`(q−1)q^{−1}` は整数でない）を明示させる。

---

## ★★Λ7 の節点表を実態に合わせて作り直す（2026-09-07。★旧表は破棄）

★★**ユーザーの指摘で発覚した**: 「Y は全部で 24 個か」という問いに答えるため数え直したところ、
★**旧表（Y1–Y16）と実際のラベルが食い違っていた。**

- 旧表の **Y15 は「LKW（Theorem 6.15）」、Y16 は「`Art(U^v)=Γ^v`」**だった
- ★**実際の Y15 は Hasse-Arf、Y16 は固定環への作用**になっている（**番号を使い回した**）
- 旧表が想定していなかった**基盤ノード**（Y14・Y16・Y17・Y20 はすべて固定環の配管）が
  実際には必要になり、そこに番号を消費した
- 枝番（`Y7a/Y7b`・`Y15/Y15b`・`Y19/Y19b+c/Y19d`）で **27 ラベル**に膨らんだ

⇒ ★★**ラベルの番号は進捗の指標として当てにならない。**
**原典の項目で数えること。** 以下を正とする。

### ★Yoshida §6 の項目（設定 2 つを除く 15 項目）

| 原典 | 状態 | 実体 |
|---|---|---|
| `def-6-1` 下付き分岐群 | ✅ | `LowerRamificationGroup`（Y1+Y2） |
| `prop-6-2` θ_0・θ_n の単射 | ✅ | `RamificationQuotientEmbedding`（Y3） |
| `cor-6-3` `e_0 ∣ n` | ✅ | `RamificationJumpDivisibility`（Y4）＋ `AbelianJumpDivisibility`（Y13、原文に忠実な別証明＋下流の受け渡し口） |
| `lemma-6-4` | ✅ | `ConjugateSumValuation`（Y5） |
| `lemma-6-5` | ✅ | `UniformizerExpansion`（Y6。★erratum 1 件目を掛け算形で回避） |
| `prop-6-6` Sen (i)(ii)(iii) | ✅ | `SenJumpFiltration`（Y7a）＋ `SenValuationCongruence`（Y7b） |
| `cor-6-7` p 進展開 | ✅ | `SenJumpExpansion`（Y8） |
| `lemma-6-8` 平均公式 | ✅ | `QuotientRamIndexAverage`（Y9） |
| `prop-6-9` **Herbrand** | ✅ | `HerbrandFunction`（Y11） |
| `lemma-6-10` 閉じた形と合成則 | ✅ | `HerbrandComposition`（Y12） |
| ★**`thm-6-11` Hasse-Arf** | ✅ **仮定ゼロ** | `HasseArf`（Y15）→ `HasseArfInduction`（Y15b）→ `HasseArfStrongInduction`（Y21）→ `SubgroupActionBridge`（Y22）→ `TameQuotientTower`（Y23） |
| `def-6-12` 上付き番号付け | ✅ | `UpperRamificationGroup`（Y18） |
| `cor-6-13` (i)(ii) | ✅ | 同上 |
| `cor-6-13` (iii) | **走行中** | Y24 |
| ★**`prop-6-14`** | **未着手** | ★Λ6 §4（`prop-4-4`(ii)(iii)・`lemma-4-3-ii`）待ち |
| ★★**`thm-6-15` LKW** | **未着手** | ★**Λ7 の目的地** |

⇒ ★**15 項目中 13 が着地、1 が走行中、残り 2。**

### ★Yoshida の項目に対応しない基盤ノード（旧表に無かったもの）

`FixedRingRamificationIndex`（Y10、`e(K′/K″)=|H|`）/ `FixedRingTower`（Y14）/
`FixedRingAction`（Y16）/ `FixedRingMonogenic`（Y17）/ `FixedRingBaseAlgebra`（Y20）——
★**この 5 本はすべて「固定環の配管」**で、Λ7 の鎖を原典の設定だけに載せるために要った。
`AbsGalRamificationFiltration`（Y19）/ `InertiaReduction`（Y19b+c）/
`RamificationFiltrationBuild`（Y19d、走行中）—— ★**pGC §2 の入力を作るためのもの**で、
これも旧表に無い。

### ☆番号の振り直しについて

★**ユーザーに判断を委ねた**（勝手に振り直すと過去の記録との対応が切れるため）。
★**当面は「原典の項目で数える」を正とし、Y ラベルは持ち場の識別子としてのみ使う。**

### Y19d 着地（部分）—— ★「完全分岐の穴」を抜けた。`compat` 以外は全部埋まった

`Found/PGC/RamificationFiltrationBuild.lean` **670 行、sorry 0**（見積 400–900 の内）。
`lake build` 成功(3199 jobs, 9.7 秒)。★`lean_start` **0 回**（12 回連続）——
共有 REPL の `lean_check` 14 回 + `leanfile.mjs` 13 往復。
`#print axioms` すべて `[propext, Classical.choice, Quot.sound]`。

### ★★★Y19 が「代入できていない」とした障害が消えた

Y19 は「Y18 の**仮定 10 個**が各有限段で完全分岐を要求するため、
橋（`comap_compat_of_coe_mul_coe_eq`）だけ先に置いた」と書いていた。
★★**Y19d が 10 個すべてを供給した**（証拠は `stage_upperRamification_coe_mul_coe_eq`）。

★★**穴の出口は Teichmüller の抽象核 1 本だった**:
`exists_fixed_sub_mem_maximalIdeal`（「`H` が剰余体に自明に作用するなら
`B^H` の剰余体＝`B` の剰余体」、**分岐・付値・Galois の語彙ゼロ、0.12 秒・一発**）。
そこから `exists_sub_mem_maximalIdeal_inertiaFixedRing`（★`L/L₀` が完全分岐）と
★★★`exists_uniformizer_adjoin_inertiaFixedRing_eq_top`（★**`𝒪_L = 𝒪_{L₀}[π]`**）が出た。
★実装者いわく **`hadj` の供給が要だった**。

### ★★近似で埋めなかった（規律が守られた）

★**持ち場に「近似で埋めて『できた』と言わないこと」と明記した結果、
Y19d は「`compat` は埋まらなかった — 近似では埋めていない」と正直に報告した。**
★さらに **`exists_absGalStage_ne_inertiaStageFiltration`** を作り、
★**自分が作ったものが近似 witness と値が食い違うことを証明した**
（`v ≤ 0` で `I_K ⊔ N`、`v` 大で `N` ⇒ `inertiaStageFiltration` とも
`trivialStageFiltration` とも異なる）。
☆**「本物である」を主張するのではなく「近似ではない」を証明する形にしたのが良い。**

### ★本体の見当の当否

- 見当 1（`lowerRamificationGroup_inertiaGal_zero_eq_top` が §6.1 の設定を与える）→ **当たり**。
  ★ただし**実際に要ったのは `G_0(L/L₀) = ⊤` ではなく、
  `Gal(L/L₀) ≤ G_0` から Teichmüller で作る `hres`** の方だった。
- 見当 3（Y19 の `comap_compat_of_coe_mul_coe_eq` が橋）→
  ★**橋としては正しいが渡れない**。橋の手前で「不分岐底変換での上付き番号付けの不変性」が要る。
- Y19b+c が残した穴のうち **逆包含 `G_0 ≤ Gal(L/L₀)` は要らなかった**
  （使ったのは `Gal(L/L₀) ≤ G_0` の側だけ）。
  ★**合成体の Galois 群の同型は要る。これが `compat` を止めている。**

`lean-idioms.md` **#134**（名前 4 連発）・**#135**（構造体の引数に依存する項を
`rw [← h]` すると motive が壊れる）・**#136**（`dite` の条件が `Nonempty`）・
★**#137**（**Y18 の右辺 `upperRamificationGroup G ϖ m` は `C` の上で計算した群**）。

☆逸脱 1・2 が誠実: 一意化元 `π` と原始元 `x` を**選択で 1 つ選んでおり、
非依存性は未証明**と明記している。

### 新しく配った持ち場: Y19e `compat` の心臓

★Y19d が名指しした残り 2 本:
1. `C = (𝒪_{L′})^{Gal(L′/L·L′₀)} ≅ 𝒪_{L·L′₀}`（Galois 降下）
2. ★**不分岐底変換で上付き番号付けが不変**（環同型＋群同型に沿った輸送）

★★**この 2 本が入れば `ramificationFiltrationOfCompat` に代入するだけで
`RamificationFiltration p` の本物が立ち、`Skeleton/PGC/Section2.lean` の
`prop_2_1` / `prop_2_2` が動く。**

★**mathlib を先に引かせる**（`Algebra.IsInvariant` / `FixedPoints` 系。
★`Mathlib.RingTheory.Invariant.Basic` は本木のどこからも import されていないと Y16 が実測しており、
**import 1 行で使えるかもしれない**）。
★**近似で埋めない規律**と、**Y19d の `exists_absGalStage_ne_...` に倣って
「近似ではない」を証明すること**を明記した。

### Y24 着地 —— Corollary 6.13 (iii)。★本体の見当が 4 つとも当たった（今日初）

`Found/PGC/UpperRamificationIndex.lean` **395 行（うち証明は 90 行弱）、sorry 0**。
`#print axioms` 11 宣言すべて `[propext, Classical.choice, Quot.sound]`。
`lake build` 成功(3114 jobs, 6.5 秒、警告 0)。★`lean_start` は呼んでいない（**13 回連続**）。

主定理 `index_upperRamificationGroup_dvd`:
`(G^m).index ∣ (Nat.card (ResidueField B) − 1) * Nat.card (ResidueField B) ^ (m − 1)`。

★★**本体の見当（θ_0 で `q−1` / θ_n で `q` / 跳びは高々 `m−1` 回 / 望遠鏡積）は 4 つとも当たった。**
☆**今日は本体の見立てが 15 回訂正されているので、4 つとも当たったのは記録に値する。**
差分は 2 つだけ:
1. `|G_{i−1}/G_i|` は **`Subgroup.relIndex` で書くのが最安**
   （`Subgroup.index_eq_card` で Y3 の `thetaMulQuot`/`thetaAddQuot` の定義域とそのまま一致、
   望遠鏡積は `Subgroup.relIndex_mul_index` の帰納 **4 行**）
2. ★**跳びの数え上げに `φ_G` の「値」は不要**で、原文の `0 ≤ φ_G(i−1) ≤ φ_G(n−1) < m` は
   **`φ_G(0)=0` と `φ_G(ψ_G(m))=m` の 2 回の狭義単調性**に落ちる

★**Y23 の Hasse-Arf（仮定ゼロ版）はそのまま使えた** ——1 回呼ぶだけで、
仮定の付け足しも読み替えも不要。★**本ノードは Y23 の仮定に 1 つも足していない。**

★**抽象核が `Subgroup` と `Finset` と `ℕ` だけ**で書け（**0.03 / 0.23 秒**）、
★**具体層は `F i := G_i`, `e := q−1`, `M := m−1`, `g i := ⌊φ_G(i)⌋₊` を代入するだけ**になった。
★**#133（`letI`/`haveI` が宣言の境界を越えない）には当たらなかった** ——
抽象核が `[Group G]` 以外のインスタンスを 1 つも要求しない形に落ちたので、
★**`letI`/`haveI` を一度も書いていない。**

★**mathlib を先に引き、新しい `def` / `instance` を 1 つも書いていない。**
★★**`Nat.card_units : Nat.card αˣ = Nat.card α − 1` が `[GroupWithZero α]` だけで
有限性を要求しない**ことを見つけ、★**剰余体の有限性を仮定せずに済ませた**
（逸脱 5: 無限なら両辺 0 で自明に真）。

☆**`m = 0` の退化を明示**: ℕ の切り詰め引き算で `q^{0−1} = 1`、右辺は `q−1`、
左辺は `|G/G^0| = 1`（`index_upperRamificationGroup_zero` を別宣言で用意）。**主張は真のまま。**

☆★**`brief.mjs` の冒頭は「半分ミスリード」だった**と実装者が報告:
「既に木にある: `UpperRamificationGroup.lean`」は **(i)(ii) だけで (iii) は無い**。
★**1b が決定打**——「(iii) の 4 行が段取りそのもので、
これを読まずに書き始めていたら**跳びの回数の上界を自分で組み直していた**」（10 人目の証言）。

`lean-idioms.md` **#138**（`Subgroup.card_dvd_of_injective` は行き先の型で単一化する）・
**#139**（`Subgroup.relIndex` は `Subgroup.index_eq_card` でそのまま商の `Nat.card` になる）。

### ★Λ7 の残りは 2 項目、どちらも Λ6 §4 に塞がれている

★**15 項目中 14 が着地**（`cor-6-13`(iii) が入った）。残るのは **`prop-6-14`** と
★★**`thm-6-15`（LKW＝Λ7 の目的地）**。
構造化係が p.17 の証明を逐行で特定したところ、**Prop 6.14 が使うのは
`prop-4-4`(ii)・`prop-4-4`(iii)・`lemma-4-3-ii` の 3 つだけ**である。

⇒ ★**D28（道 A）で確定した 3 本のうち 1 本目
`[K̂^m_f : K̂^ur] = q^n − q^{n−1}` を配った。**
★在庫調査係が名指しした穴「★**`K^ur → K̂^ur` の完備化を渡る段が無い**」が
本当に中身かを測らせる（★本体は確かめていない）。
★★**Lubin-Tate 機械（57 ファイル / 14,525 行）を 1 行も触らないこと**を明記した
（D28 の「道 A は 471 宣言を 1 件も触らない」を守るため）。
★**`[Fintype (ResidueField _)]` を要求する在庫に当たったら、それが D28 の予測どおりの
詰まり所である**ことも警告した。

### ★`brief.mjs` の「既に木にある」を項目の**部分**まで断るようにした（Y24 の指摘）

★**Y24 の報告**: 「冒頭の『既に木にある: `UpperRamificationGroup.lean`』は**半分ミスリード**
（(i)(ii) だけで (iii) は無い）」。
★**照合は `data-item` の「項目」単位なので、原典の 1 項目が (i)(ii)(iii) に分かれていると
「一部だけ埋まっている」場合に誤解を招く。**

**直し方**: 逐語に**部分番号**（`(i)` `(ii)` `(iii)` `(1)` …）が **2 つ以上**見えるときだけ、
冒頭と末尾の両方に断りを足す（noise を出さないため）。

**検算**:
| 項目 | 部分番号 | 木にある | 断りが出るか |
|---|---|---|---|
| `cor-6-13` | (i)(ii)(iii) | ✅ | ★**出る**（正しい） |
| `prop-6-6` | (i)(ii)(iii) | ✅ | ★出る |
| `cor-6-7` | 無し | ✅ | ★**出ない**（正しい。noise を出さない） |

★`--json` は無変更（`selfNodes` は元から入っている）。
★**Proof 抽出に退行なし**（Yoshida08 68 項目で 31/20/0 のまま）。

☆★**道具の警告が「半分正しい」ときがいちばん危ない**、という実例。
2026-09-07 の朝は**本体が警告を切り落として**既出の項目を配り、
夕方は**警告が粗すぎて**実装者に「半分ミスリード」と言われた。
★**どちらも道具側で塞いだ**（冒頭＋末尾の二重掲示 → 部分番号の断り）。

### Λ6 §4-a 着地 —— `[K̂^m_f : K̂^ur] = q^n − q^{n−1}`（D28 の 1 本目）

`Found/PGC/LubinTateCompletionDegree.lean` **514 行（証明本体は 120 行弱）、sorry 0**。
`#print axioms` 7 宣言すべて `[propext, Classical.choice, Quot.sound]`。
`lake build` 成功(3433 jobs, 10 秒)。★`lean_start` **0 回**（14 回連続）。
往復は具体層 8・抽象核 4 の**計 12 回**、うち **6 宣言が一発**。

★**抽象核 3 本がすべて「分岐・付値・Lubin-Tate が 1 語も出ない」形に落ちた**
（`finrank_adjoin_of_isEisensteinAt` **0.12 秒・一発** /
`isEisensteinAt_map_of_maximalIdeal_eq_span` 0.26 / `adjoin_image_eq_adjoin_simple` 0.07・一発）。
**36 回連続。**

### ★★在庫調査係が名指しした穴は「半分外れていた」—— 本当の穴は別だった

在庫調査係（D28）は「★`K^ur → K̂^ur` の**完備化を渡る段が無い**」と書いた。
★**実装者の判定: 整数環のレベルでは既に木に在った**
（`maximalIdeal_unramifiedCompletionInt_eq_span` + `baseIntHom_eq_uniformizerCompletionInt`）。
足りなかったのは**その等式を Eisenstein 多項式に載せる 1 本**（抽象核、**20 行**）だけで、
★**局所準同型も DVR も要求しない**（単元が単元へ写るだけで足りる）。

★★**本当の穴は `Λ_n ⊆ K(α)`（原典 (ii) の `µ_{f,m} ⊂ L(α)`）で、
`L = K` の場合ですら木に無かった。** これが最大の未知だったが、
`iteratedLubinTateTorsionPoints_eq_union` + `iteratedLubinTatePsiTorsionPoints_subset_adjoin` +
`lubinTateActionAtTorsionPoint_pi_mem_…` の `n` 帰納法 **40 行**で閉じた
（★決定打は **`lubinTateActionAtTorsionPoint` の値がそもそも `adjoinIntegers K α ⊆ K(α)` にある**こと）。
★★**原典 (i)（`𝒪/𝔭^m ≅ µ_{f,m}`）を作らずに済んだ。**

### ★D28 の予測は当たった

- ★**`[Fintype (ResidueField _)]` で詰まらなかった** —— `𝒪_{K̂ur}` について `Fintype` を
  要求する在庫を**一度も使っていない**（`isDiscreteValuationRing_` と
  `maximalIdeal_…_eq_span` の 2 本だけで足りた）。
- ★★**Lubin-Tate 機械 471 宣言を 1 行も触っていない**（道 A の約束が守られた）。
- ★**`L = K` 版の次数は流用できず書き直しになった** ——
  あちらは `Gal(K(α)/K) ≃ (𝒪/π^n)^×` から次数を読む道で、底を替えると
  **Galois 群の同定からやり直し**になる。本ノードは**多項式の既約性から直接**読んだ
  （Eisenstein → Gauss → `minpoly`）。★**原典 p.7 の証明（付値の挟み撃ち）より短い。**

★**1b が決定打**（12 人目の証言）:
> ★**1b を読んで初めて「原典は付値の挟み撃ちで既約性を出している」と分かり、
> 木の Eisenstein 在庫で置き換えれば付値を一切通らずに済むと判断できた。
> 1b が無ければ原典の段取りを組みに行って数十往復していた。**

`lean-idioms.md` **#140**（`Polynomial.Monic` の項に**ドット記法は使えない**）・
★**#141**（★**`lean/ABC3/Found.lean` は CRLF**。Python の置換は `'rb'`/`'wb'` で。
`python - <<'PYEOF'` も PreToolUse フックに潰される＝#117(iv) の Python 版）。

### ★★★本体の持ち場の誤り（今日 16 件目）—— `pdfPage` を手で書いて間違えた

本体は持ち場に `pdfPage := 8` と書いた。★**誤り。正しくは 7。**
検算: 構造化 HTML の `data-pdf-page="7"` / `brief.mjs` の雛形も **7**。
★**実装者が「指示が誤りだと思われる」と申告し、機械抽出値の 7 を採った。正しい判断。**

☆★**原因は本体が雛形を写さず手で書いたこと。**
★**規則: `.src` は `brief.mjs` の雛形をそのまま写す。手で書き直さない。**
（M31 が「G1 は `.src` の `pdfPage` と HTML の `data-pdf-page` の一致を見ていない」と
測っているので、★**この種の誤りはゲートで捕まらない。**）

### 新しく配った持ち場: Λ6 §4-b `Gal(K̂^m_f/K̂^ur) ≃* (𝒪_K/π^n)^×`（D28 の 2 本目）

★Λ6 §4-a の実装者が「**本ノードの `lubinTateCompletionField_eq_adjoin_simple` と
`finrank_lubinTateCompletionField` をそのまま入力にできる**」と書いているので、そのまま配った。
★**Λ6 §4-a が「底を替えると Galois 群の同定からやり直しになる」と指摘した、その同定が中身。**
★`.src` は**雛形をそのまま写すこと**を明記（本体の誤りの再発防止）。
★**#69（`adjoinField`/`adjoinIntegers` の境界は 212 秒で timeout）に当たる公算が高い**ことを警告
（`Λ` の作用＝整数側と `Gal`＝体側を行き来するため）。

### ★★測定: `.src` の `pdfPage` を機械で検査できるか —— **今日の誤りは捕まらない**（負の結果）

本体が `pdfPage := 8` と書いて誤った（正しくは 7）ので、
★**同じ種類の誤りをゲートで捕まえられるか**を測った（★`check.mjs` は書き換えていない）。

**前提**: M32（メタ第 11 回）は「`.src` の `pdfPage` と HTML の `data-pdf-page` の**完全一致**」を
測って **誤報 985 件**と判定し、★**頁の検査を意図的に外した**
（項目が数頁にまたがるとき、`.src` は引いている段の頁を、HTML は項目の先頭頁を指すため）。

**本測定が試したのはより弱い条件**: 「`.src` の `pdfPage` が
**その項目の頁の範囲**に入っているか」（範囲 = [項目の `data-pdf-page`,
同じファイルで次に頁が進む項目の `data-pdf-page`]）。

**結果（`.src` 総数 4,644、うち構造化に id があるもの 4,641）**:

| ずれ | 件数 |
|---|---|
| −3 | 1 |
| −2 | 2 |
| **−1** | **56**（★系統的。見出しが前頁の末尾から始まる） |
| **0（範囲内）** | **4,576** |
| +1 | 2 |
| +3 | 4 |

⇒ ★**規則を「[lo−1, hi] に入ること」にすると誤報は 6 件**（重複除去後）で、
M32 の 985 件より **2 桁小さい**:
```
FrdI#frdi-def-2-4        pdfPage=51 範囲=[47..48] ずれ +3  （3 ファイル）
GenEll#genell-lemma-3-5  pdfPage=15 範囲=[17..18] ずれ −2
GenEll#genell-prop-1-4   pdfPage=3  範囲=[6..8]   ずれ −3
GenEll#genell-prop-1-4   pdfPage=9  範囲=[6..8]   ずれ +1
```

### ★★★しかし本測定の結論は「この検査では今日の誤りは捕まらない」

★**`Yoshida08#prop-4-4` の範囲は `[7, 8]` と計算される**（次に頁が進む項目が p.8 にあるため）。
⇒ ★★**本体が書いた `pdfPage := 8` は「範囲内」と判定され、捕まらない。**

☆★**負の結果を正直に記録する。** 検査を作っても動機となった誤りは防げない。
★**頁の 1 つのずれを機械で捕まえるには PDF 本文まで読む必要がある**が、
`check.mjs` は既に逐語照合で `pdftotext` を叩いており、そこまでやると
**ゲートが重くなる**（キャッシュがあっても 7 秒 → 55 秒の実績）。

⇒ ★**本当の対策は道具ではなく手順**である:
★★**`.src` は `brief.mjs` の雛形をそのまま写す。手で書き直さない。**
（本体の誤りは「雛形を見ずに手で書いた」ことが原因だった。実装者は雛形を採って正した。）

☆**副産物**: 上の 6 件は**実際に locator がずれている可能性がある**。
★**急ぎではない**（G1 の項目名照合は通っているので、頁だけの問題）が、
**人が見る価値のある候補**として記録する。

☆**再測定のコマンド**（本体が再実行できる形）:
`node <scratchpad>/outliers.mjs`（測定スクリプトは scratchpad。
★**ヒアドキュメントでは書けない** —— バックスラッシュが食われる（#117(iv)）。
**Write ツールで書くこと**。本測定で実際に 1 度踏んだ）。

### Λ6 §4-b 着地 —— `Gal(K̂^m_f/K̂^ur) ≃* (𝒪_K/π^n)^×`（D28 の 2 本目）

`Found/PGC/LubinTateCompletionGalois.lean` **635 行・宣言 22 本、sorry 0**。
`#print axioms` 5 本すべて `[propext, Classical.choice, Quot.sound]`。
`lake build` 成功(3436 jobs、自モジュール 16 秒)。★`lean_start` は呼んでいない（**16 回連続**）。
抽象核 7 本は `lean_check` 4 往復（★**4 本まとめて 0.30 秒・一発**）、具体層 15 本は
`leanfile.mjs` **2 往復**。

### ★★★本体の見立てが 2 つとも外れた（今日 17 件目・18 件目）—— どちらも良い方向

1. **「`L = K` 版の骨格を `L = K̂` で書き直す仕事」→ ★丸ごと流用できた。**
   ★**原典の 2 段（`Gal ≅ Aut_𝒪(µ) ≅ (𝒪/𝔭^m)^×`）を経由しない**設計にし、
   `galoisCompletionEquivBase : Gal(K̂^ur(α)/K̂^ur) ≃* Gal(K(α)/K)` を作って
   **既存 425 行の `galoisReciprocityEquiv` をそのまま合成**した。
   ★★**安くなった理由: `minpoly_{K̂^ur}(ι α) = ψ_n = minpoly_K(α)`
   （Eisenstein が完備化を渡るので**底を替えても最小多項式が変わらない**）。**
   ⇒ 単射は「`x` での値で決まる」だけ、全射は `PowerBasis.equivOfMinpoly` だけ。
   ★本体が挙げた在庫のうち `unitActionQuotientLift` / `lubinTateActionAtTorsionPoint` /
   `principalUnitsQuotientEquiv` / `reciprocityMap` は **1 つも直接使っていない**。
2. **「#69（`adjoinField`/`adjoinIntegers` の境界）に確実に当たる」→ ★当たらなかった。**
   ★**整数側に一度も降りない設計**にしたため（`Λ` の作用は `galoisReciprocityEquiv` の
   内側に閉じ込め、体側の在庫だけで回した）。

★**Λ6 §4-a の出力はそのまま入力にできた**（`lubinTateCompletionField_eq_adjoin_simple` を 2 回）。
★**`finrank_lubinTateCompletionField` は使わなかった**（次数を経由せずに済んだ）。
★**`[Fintype (ResidueField _)]` に詰まらなかった**（§4-a と同じく一度も使わずに済んだ）。
★**mathlib を先に引き、「無い」と判定した項目はゼロ**（5 波連続）。

☆**`brief.mjs` の冒頭が効いた**（今日入れた部分番号の断りが初めて役に立った）:
> 冒頭の「既に木にある／部分番号 (i)(ii)(iii) に注意」で
> ★**(ii) の次数だけが埋まっていて (iii) は空**と即断できた。
★**1b も効いた**: 「原典も単射＋濃度で押している」と分かり、
**こちらは濃度の代わりに最小多項式で全射を直接出す**（濃度計算も Galois 性も不要）判断ができた。

☆逸脱 3 が誠実: `galoisCompletionReciprocityEquiv` は**原始点 `α` の選択に依存する**
（原典の `ρ_{f,m}` は `α` に依らない）。★**選択を隠した形 `nonempty_…` を別に用意した。**

`lean-idioms.md` **#142**（`set` は他の仮説の型に現れる項を抽象化するとその仮説を `τ✝` に化けさせる。
1.72 秒 → 0.26 秒）・**#143**（`∃!` に `refine ⟨_, ?_, ?_⟩` すると β 簡約されず `rw` が落ちる。
★**一意性側は `intro` が β 簡約するので非対称**）・**#144**（複数行の `calc` の第 1 項が
関数適用だとパーサが途中で切る）。

### ★測るべき重複（§4-b の実装者が指摘）

`galoisCompletionEquivBase`（`Gal(K̂^ur(α)/K̂^ur) ≃* Gal(K(α)/K)`）は
★**汎用の「不分岐底変換」**で、★**同時に走っている `UnramifiedBaseChangeInvariance`（Y19e）と
主張が近い可能性がある。** ⇒ 次の持ち場（§4-c）に「重複を測ること。走行中なら
『測れなかった』でよい。★**待たないこと**」と書いた。

### 新しく配った持ち場: Λ6 §4-c `N(−α) = π` と「α が素元」（D28 の 3 本目）

★**D28 の 3 本のうち 2 本が今日着地した**ので、残り 1 本。
★あわせて Prop 4.4(ii) の残り「α が素元」も含めた（原典が同じ 1 文で述べている）。
★★**在庫調査係の測定: `Algebra.norm` を型に含む宣言は木全体で 4 件、`Found/PGC` には 0 件。**
⇒ ★**ノルムは本当に新しい領域である。** mathlib を先に引かせる。
★本体の見当「`N(−α)` は最小多項式の定数項に落ちる」を
**「実機で確かめていない。外れたら訂正して報告すること」**と明記して配った
（★**本体の見立ては今日 18 回訂正されている**と数字も添えた）。

---

## ★★★★★`RamificationFiltration p` の**本物が仮定ゼロで立った**（Y19e、2026-09-07）

`Found/PGC/UnramifiedBaseChangeInvariance.lean` **960 行 / 54 宣言、sorry 0**。
`#print axioms` **14 宣言すべて `[propext, Classical.choice, Quot.sound]`**。
`lake build` 成功(3212 jobs, 17 秒)。★`lean_start` **0 回**（共有 REPL が使用中だったため
`leanfile.mjs` のみ、約 27 往復）。

```lean
noncomputable def ramificationFiltration (p : ℕ) [Fact p.Prime] : RamificationFiltration p
-- ★★仮定は 1 つも無い
```

★**`Interface/PGC/LocalFieldData.lean:160` の構造体は
`Skeleton/PGC/Section2.lean` の `prop_2_1` / `prop_2_2` が入力に要求しているもの**で、
`waiting` に「**Herbrand の定理を要する**」と書かれていた。
★★**Herbrand は今日 Y11 が着地させ、5 波（Y18 → Y19 → Y19b+c → Y19d → Y19e）で
本物が立った。**

### ★★退化していないことを 3 通りで示した

木には退化 witness が 2 つある（Y19 の `trivialStageFiltration`＝下からの `⊤`、
Y19b+c の `inertiaStageFiltration`＝上からの `I_K`）。
★**Y19e は自分の作ったものが両方と値が食い違うことを証明した**
（`exists_stage_ne_inertiaStageFiltration` / `stage_ne_trivialStageFiltration`）。
★★さらに **`coe_ramificationFiltration_mul_coe`**（`Γ_K^v · M = Gal(L/L₀)^v` の引き戻し。
**退化した族では成り立たない式**）を作った。
☆**「本物である」を主張するのではなく「近似ではない」を 3 通りで証明する形**が守られた。

### ★★★本体の見当が 3 つとも訂正された（今日 19・20・21 件目）

1. 「環同型＋群同型に沿った輸送」→ ★**半分外れ**。要るのは
   ★**全射準同型に沿った輸送**（`I(L′/K) ↠ I(L/K)` は**同型ではない**）。
   ★**`ramIndex` さえ対応すれば `G_n`・`φ`・`ψ`・`G^m` が自動で対応する**、という形に落ちた。
2. 「(1) 固定環＝底の整数環が要る」→ ★★**要らなかった。**
   `C = (𝒪_{L′})^H` を `𝒪_{L·L′₀}` と同定せず、
   **`𝒪_L → C` の埋め込みだけ作って「素元が素元のまま」を示す道**に変えたので
   **(1) は 1 行も使っていない**（★段取りとの最大の差分）。
3. ★★**本体が予期していなかった真の残りは `e = |I|`** だった
   （`ramificationIndex_eq_card_inertiaGal`）。これも埋めた。決め手は
   `AbelianSplitUnramified.lean` の `exists_unramified_subextension`（最大不分岐部分拡大の次数が `f`）
   ＋ Y19b+c の `lift_fixedField_inertiaGal`。

★**mathlib を先に引き、11 本を採用**。★**`Nat.sSup_le` と `IntermediateField.finrank_lift` は
測って不在**（`#check` で `Unknown constant`）と正しく判定し、`csSup_Iic` と
`IntermediateField.equivMap` + `LinearEquiv.finrank_eq` で代替した。
★`Algebra.IsInvariant` は**引いたが使わなかった**（(1) が不要になったため）。

`lean-idioms.md` **#145**（`adjoinIntegersIncl` の係数を 2 層いっぺんに `rfl` で潰すと
**kernel deterministic timeout**）・
★★**#146**（**`ker_restrictNormalHom_eq_fixingSubgroup` は同名が 2 つある** ——
`AbsGalRamificationFiltration.lean:492` と `LubinTateClosure.lean:127`。
両方が import 圏に入ると曖昧で落ちる → mathlib の `IntermediateField.restrictNormalHom_ker`）・
**#147**（mathlib の `Subgroup.mul_normal` は `↑(H ⊔ N) = ↑H * ↑N` の向き。教科書の記憶と逆）。
☆`set_option maxHeartbeats` を 3 宣言で使用（逸脱 6 に明記）。

### ★★`prop_2_2` の 3 要求のうち 2 つが今日埋まった

`Skeleton/PGC/Section2.lean` の `prop_2_2` が挙げる「依拠する境界外の結果」:

| 要求 | 状態 |
|---|---|
| `RamificationFiltration p` | ★★**Y19e が今日作った（仮定ゼロ）** |
| 上付き↔下付き番号付けの変換（Serre Ch. IV） | ★**Y18 が作った** |
| ★**`Γ_K^0 = I_K`**（Cor 1.3 の系） | **残り。Y19f として配った** |

☆★**`prop_2_1` も近い**: 3 つの要求のうち
**p 進対数は解消済み**（`padicLog_bijOn`、sorry 無し）、
**Verlagerung は mathlib にある**（`MonoidHom.transfer`）。
残るのは `.implicitStep`「log による `U_K⊗Q_p ≅ K` と Verlagerung の両立性から
Prop 2.1 自体への一段」だけ。

### 新しく配った持ち場: Y19f `Γ_K^0 = I_K`

★**入口は Y19d の `absGalStage_of_nonpos`（`v ≤ 0` で `I_K ⊔ N`）**のはず、と見当を書いた
（★**実機で確かめていない**と明記）。
★抽象核は「開正規部分群の基本近傍系を走る `⋂_N (H ⊔ N) = H`（`H` が閉のとき）」という
**純位相群論**に落ちるはず。
★★**退化 witness で「出した」ことにしないこと**を強く書いた ——
`inertiaStageFiltration` なら `Γ_K^0 = I_K` は**自明に成り立ってしまう**ため。
★**#146（同名 2 つ）に確実に当たる**ことも警告した（両方を import する圏に入るので）。

### Y19f 着地 —— `Γ_K^0 = I_K`。★★本体の見当が当たり、18 回続いた訂正が止まった

`Found/PGC/RamificationFiltrationZero.lean` **356 行 / 16 宣言、sorry 0**。
`lake build` 成功(3224 jobs, 6.4 秒)。★`lean_start` **0 回**（18 回連続）。
`leanfile.mjs` **3 往復**、`set_option maxHeartbeats` は**不要だった**。

主定理 `ramificationFiltration_Gv_zero : (ramificationFiltration p).Gv K 0 = absInertia K`。

★**抽象核が `ContinuousMul` だけで済んだ**（★`N` の開性すら要らない）:
`Subgroup.iInf_sup_eq_self_of_isClosed`（**0.19 秒**）。
★さらに**一般形** `Subgroup.iInf_sup_eq_topologicalClosure`（`H` に**何も仮定せず**
`= H.topologicalClosure`）を **0.14 秒・一発**で出し、
★**必要性** `Subgroup.isClosed_iInf_sup`（左辺は常に閉 ⇒ `IsClosed H` は落とせない）まで作った。
☆**「仮定が要る」ことを自分で証明する形**は良い作法。

★**本体の見当「`absGalStage_of_nonpos` を逆極限に通せば出る」はそのまま当たった。**
☆**今日 18 回連続で続いていた本体の訂正が、ここで止まった。**

### ★退化 witness で「出した」ことにしない規律が 3 通りで守られた

1. 主定理の主語が `(ramificationFiltration p).Gv K 0` そのもので、
   ★**証明が `ramificationFiltration` の定義を経由することを型で強制**している。
2. ★**`trivialStageFiltration_limit_zero_ne_absInertia`** —— 退化 witness では
   `Γ^0 = ⊤ ≠ I_K` で**本定理は偽**。
   ★**`absInertia_ne_top` を無条件で証明**（`Γ_K/I_K` が任意の `ℤ/n` へ全射することから）。
   ⇒ ★**主定理は「どんな段データでも成り立つ空虚な主張」ではない。**
3. `exists_coe_ramificationFiltration_mul_coe_ne_inertiaStage` ——
   ★**`inertiaStageFiltration`（`Γ^0 = I_K` が自明に成り立つ上からの近似）とも値が食い違う。**

### ★★本体の持ち場に対する訂正 1 件（実装者の測定）

本体は「`Skeleton/PGC/Section1Cor13.lean` の `inertia_recoverable` から `Γ_K^0 = I_K` が出るか」
と書いた。★**測った結果 NO。両者は論理的に独立。**
`inertia_recoverable` は「`I_K` が `Γ_K` から**群論的に復元できる**」であって
★**`Γ_K^0` の値については何も言わない。**
`Γ_K^0 = I_K` には分岐論の入力（`v ≤ 0` で上付き分岐群＝惰性群）が要る。
★**Corollary 1.3 が効くのは後半**（「`Γ_K^0 = I_K` なので `Γ_K^0` が復元できる」の側）で、
それを `ramificationFiltration_Gv_zero_recoverable` として別宣言に置いた。
☆★**1b が訂正の入口だった**: 「原典はこの項目に Proof を付けていない/論証が主張の**手前**にある」
で cor-1-3 の地の文が「不分岐判定 `q_L = q^{[Γ_K:H]}`」であって `Γ_K^0` の値ではないと分かった。

### ★★`prop_2_2` の名指しの 3 入力が全部そろった

| 要求 | 実体 |
|---|---|
| `RamificationFiltration p` | `ramificationFiltration`（Y19e、**仮定ゼロ**） |
| 上付き↔下付き変換 | `upperRamificationGroup_eq_lowerRamificationGroup`（Y18） |
| `Γ_K^0 = I_K` | `ramificationFiltration_Gv_zero`（Y19f） |

★**ただし実装者は「`prop_2_2` 自体が閉じるとは言っていない」と正しく留保した** ——
同定理は `IntKbar` / `CompKbar` を**未構築の対象として抽象化**しており、
★**具体構成が「別途 `Found/` の課題」として残っている。**

☆`#146`（同名 2 つ）に**当たらなかった**理由も報告された:
`ker_restrictNormalHom_eq_fixingSubgroup` を一度も書かず、木の包み
`ker_restrictNormalHom_unramifiedClosure` を引いたため。
★**import 圏には両方入っていたので、名前を書いていたら確実に落ちていた。**
`#147` には当たったが、★**先に idiom を読んでいたので往復 0 で済んだ。**
`lean-idioms.md` **#148**（`∀ U ∈ s, P U` の穴に `{U}` を暗黙にした補題は嵌まらない）。

### 新しく配った持ち場: Y25 `𝒪_{K̄}` と `K̄^∧` を `Γ_K`-加群として構成する

★★**`prop_2_2` の docstring の前提が今日の M3 で古くなっていた**（本体が実測）:
docstring は「スペクトルノルムは**有限次拡大にのみ適用**なのでそのままでは使えない」と書くが、
★**`closureNormedField`（`K.closure` 全体のノルム）と
`closureCompletion K`（= `K̄^∧`）と `closureCompletionInt` が既にある**
（★**M3 は今日着地したばかり**）。`absGal` と `closure` を同時に含む宣言は **65 件**。

⇒ ★**`CompKbar K := closureCompletion K` はそのまま使えるはず**、
★**`IntKbar K := {x : K.closure // ‖x‖ ≤ 1}`**（★`closureCompletionInt` は
**完備化の**整数環なので別物。`K̄ ⊊ ℂ_K`）と見当を書いた（★**実機で確かめていない**と明記）。
★★**自明な作用（`σ • x = x`）で型クラスを満たしたことにしないこと**を強く書いた ——
`DistribMulAction` は自明作用でも満たせてしまうため。

### Λ6 §4-c 着地 —— ★★★原典の等号が木の設定では**偽**だと分かった（4 件目の食い違い）

`Found/PGC/LubinTateCompletionNorm.lean` **444 行、sorry 0**。
`lake build` 成功(3437 jobs、自モジュール 15 秒)。
★`lean_start` **1 回だけ**（17.3 秒。見積 90 秒より速い）、以後 `lean_check` **11 回**、
`leanfile.mjs` は不要だった。

### ★★★`N(−α) = π^{ϕ^{m−1}}` は字面のまま偽 —— 原典と木で `f` の型が違う

★**原典は `f ∈ 𝒪_L[X]`（多項式）に限っている**（逐語 "Let m ≥ 1 and f ∈ O_L[X] as above"、
Lemma 4.3(i) の証明中 "as f_m is a monic in O_L[X]"）。
そのとき `f_m/f_{m−1}` は**多項式の商**で定数項がちょうど `π^{ϕ^{m−1}}`。

★★**ところが木の Lubin-Tate 機械は `f : PowerSeries 𝒪_K` を扱い、
`ψ_n` は `[π^n]/[π^{n−1}]` の Weierstrass distinguished 部分である。**
木の `iteratedLubinTatePsi_coeff_zero_mul` は `ψ_n.coeff 0 · U'_n(0) = π` としか言っておらず、
★**木自身も `norm_iteratedLubinTatePsi_coeff_zero` で「ノルムが等しい」までしか主張していない。**

★**実装者が反例を挙げた**（手計算、Lean には書いていないと明記）:
`K = ℚ_2, π = 2, q = 2, f = 2X + X² + 2X³`（★木の 3 仮定をすべて満たす）。
`r_1 = 2 + X + 2X²`、`𝔪` 内の根は `α ≡ 22 (mod 32)`、`ψ_1 = X − α`、
★**`ψ_1(0) = −α ≡ 10 (mod 32) ≠ 2`**。`[K̂^1_f : K̂^ur] = q−1 = 1` なので `N(−β) = −α ≠ π`。
★**`v(−α) = 1`（素元性）は保たれ、単元倍だけずれる。**

⇒ ★**3 段で述べた**（statement をねじ曲げず、正しい形に直した）:
(1) `N(−β) = ψ_n.coeff 0`（**無仮定・厳密**）/
(2) `N(−β) = π·u`（`u ∈ 𝒪_K^×`）・`‖N(−β)‖ = ‖π‖`・`(N(−β)) = (π)` /
(3) `U'_n(0) = 1` を仮定して**原典の字面** `N(−β) = π`。

☆★**この論文で見つけた食い違いは 4 件目**だが、性質が違う:
1・2 件目は**誤植**（Lemma 6.5 の 1 ずれ / Prop 6.6 の "Let j = 1"）、
3 件目は**論理の穴**（Thm 6.11 の帰納法の仮定が引けない、Y21 が発見）、
★**4 件目は「原典と我々で設定が違う」**（原典は多項式、木は冪級数）。
★**原典に誤りは無い。我々の方が一般な設定を採っている。**

### ★本体の見当の当否

- 「`N(−α)` は最小多項式の定数項に落ちる」→ ★**当たった**
  （`PowerBasis.norm_gen_eq_coeff_zero_minpoly` 系でそのまま）
- 「定数項が `π`」→ ★**外れた**（`π` の単元倍まで）
- 「`π^{ϕ^{m−1}}` は `L = K̂^ur` で `π` に落ちる」→ ★**当たった**（`π ∈ 𝒪_K` は `ϕ` で固定）
- ★**「α が素元」は閉じた**（`‖α‖^[K̂^m_f:K̂^ur] = ‖π‖`）

★**抽象核の中身が「原典が `N(α)` でなく `N(−α)` と書く理由」だった**:
「`N(−1) = (−1)^{dim}` と mathlib の `(−1)^{dim}·coeff 0` が打ち消し合う」だけ（**0.18–0.33 秒**）。
☆**原典の記法の理由が形式化で説明された例。**

★**§4-b と Y19e の重複は「無し」と測れた**（`grep -c "unramifiedCompletion\|closureCompletion"
`UnramifiedBaseChangeInvariance.lean` → **0**）。Y19e の `restrictGalHom` は**同じ底 `K` 上での
`L ≤ L'` の制限**、§4-b は**底の取り替え `K → K̂^ur`**。★**別物。**

★**#69 に当たらなかった** —— 中間体の整数環 `𝒪_{K̂^m_f}` を一度も作らず、
「素元」を**ノルムの言葉**で書いた（§4-b と同じ「整数側に降りない設計」が 2 波続いた）。
★**mathlib を先に引き、「無い」と判定した項目は 0**（8 波連続）。
★本体の在庫表に無かった `norm_iteratedLubinTatePsi_coeff_zero` /
`spectralNorm_root_iteratedLubinTatePsi`（`LubinTatePsiNorm.lean`）が
**「素元」を 2 往復で片付けた**。

★**1b が今回の最大の当たり**（実装者の証言、13 人目）:
> 「the constant term of `f_m/f_{m−1}` reads `π^{ϕ^{m−1}}`」を読んで初めて
> ★**原典が多項式の商を取っていると分かり、木の `ψ_n` とのズレに気づけた。**
> ★**1b が無ければ `ψ_n.coeff 0 = π` を証明しようとして無限に往復していた。**

`lean-idioms.md` **#149**（`π` は `hπmax` の型に現れるので `rw [← h]` は必ず
`motive is not type correct`。移項補題に切り替えると 1 往復）。

### ★D28（道 A）の 3 本が全部閉じた

| # | 内容 | 実体 |
|---|---|---|
| 1 | `[K̂^m_f : K̂^ur] = q^n − q^{n−1}` | `LubinTateCompletionDegree.lean`（514 行） |
| 2 | `Gal(K̂^m_f/K̂^ur) ≃* (𝒪_K/π^n)^×` | `LubinTateCompletionGalois.lean`（635 行） |
| 3 | `N(−α) = π`（正しい形）・α が素元 | `LubinTateCompletionNorm.lean`（444 行） |

★★**Lubin-Tate 機械 57 ファイル / 14,525 行を 1 行も触っていない**（3 波とも D28 の約束を守った）。

☆新ノード候補（★**Λ7 の閉路には入っていない**と実装者が判定）:
「`f` が多項式なら `U'_n(0) = 1`」——これが埋まれば原典の字面が完全に復元する。
★**数学が足りないのではなく、木が原典より一般の `f` を扱っているために生じた節点。**

### ★★Prop 6.14 には `Definition 5.3`（§5）も要ることが分かった

`prop-6-14` の逐語: 「Let `L = K_n` and `K^m_x` **as in Definition 5.3**」。
⇒ ★**Λ7 の終盤は §4 だけでなく §5 の定義も要る。** 本体は当初「§4 の 3 つだけ」と見ていた
（構造化係の測定は「Prop 6.14 が使う **§4 側の**依存先は 3 つだけ」であって、
★**§5 側は数えていなかった**）。★**新規配布を止め、ゲートを優先する。**

---

## ★★★★ゲート前検査で **NG 42 → 13** —— 29 件は道具の取りこぼしだった（2026-09-07）

`Y25` が走行中でも `check.mjs` は **lake を使わない**ので、ゲート本体（全体ビルド）の前に回した。

### 経過

| 段階 | NG | 内容 |
|---|---|---|
| 初回 | **44** | |
| 本体が 2 件直した後 | **42** | 今日の Lean 由来の NG は 0 になった |
| ★**実体表を広げた後** | ★★**13** | ★**29 件は偽の NG だった** |

### ★本体が直した 2 件（今日の Lean 由来。実体は「逐語の写し方」）

1. `SubgroupActionBridge.lean:331` —— `Z[bb]_{≥0}` と書いていたが
   ★**正しい docstring 用の形は `Z[bb]_≥0`（波括弧なし）**。
   `brief.mjs` の「docstring に貼る逐語」と 1 文字ずつ突き合わせて直した。
2. `UpperRamificationIndex.lean:298` —— `G_{i−1} ≠ G_i` と書いていたが
   ★**`pdftotext` は `≠` の斜線を落とす**ので `G_i−1 = G_i` が正しい投影
   （[[pdftotext-drops-negation]]）。★**波括弧も出ない。**
   ⇒ 直したうえで「これは `pdftotext` の投影である」と docstring に明記した。

☆★**どちらも `brief.mjs` の雛形を写していれば起きなかった**
（★本体の `pdfPage := 8` の誤りと同じ原因）。

### ★★★29 件の偽 NG の正体 —— `check.mjs` の `ENTITIES` が標準実体を取りこぼしていた

Yoshida §2–§4 の構造化（68 件）で **S4 が 12 件落ちていた**が、
★**逐語の誤りではなく `check.mjs` が `&isin;` `&pi;` `&middot;` `&cap;` `&equiv;` `&rArr;`
などを解釈できなかった**ためだった（`decodeEntities` は表に無い実体を**そのまま残す**）。

★**これらはコーパス全体で使われている**（実測: `&middot;` 既存 55 / `&sigma;` 45 /
`&isin;` 22 / `&pi;` 29）。★**構造化係は正しく書いていたのに、道具が読めていなかった。**

☆**構造化係は正しく気づいていた**: 報告の「守れなかった・迷ったもの」に
> README への追記が §7 の「状態表」を少しはみ出した ——
> 「`.verbatim` で使える名前つき実体は `check.mjs` の `ENTITIES` にあるものだけ」という
> 失敗形の注を 1 段落足した

と書いていた。★★**つまり「道具の制約に合わせて書き方を狭める」方向で対処されており、
道具の側を直す発想が無かった。** ★**今回は道具を直した。**

**直し方**: `ENTITIES` に 30 個追加（`isin` `ni` `middot` `cap` `cup` `equiv` `ne` `empty`
`infin` `bull` `rArr` `lArr` `hArr` `Prime` `prod` `sum` `part` `radic` および
ギリシャ小文字 12・大文字 6）。

**手順（★増えたら戻す、を先に決めてから測った）**:
| 検査 | 結果 |
|---|---|
| NG 件数 | **42 → 13**（★**29 件減、増加ゼロ**） |
| `check.mjs --selftest` | ★**50/50 PASS**（器具は壊れた入力を落とせている） |
| `mojibake.mjs` | ok 文字化けなし |

★★**13 はセッション開始時の真の基準値と一致する。**
⇒ ★**今日の 26 ファイル・約 15,000 行は NG を 1 件も増やしていない。**

☆★**`&ne;` は `≠` に開くが `pdftotext` は斜線を落とす**ので、
★**`data-txt="="` を併記しないと通らない**ことをコメントに残した。

### ★残る 13 件は既存の繰り越し

G9 繰り越し **27 件**（CorrHyp の非空虚性対照が無い）と G1 繰り越し **3 件**は
**別枠で数えられている**（`NG 13` には含まれない集計）。
★**どちらも `CorrHyp/**` 由来で、D26 により本セッションでは触らない。**
