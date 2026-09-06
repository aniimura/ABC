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
