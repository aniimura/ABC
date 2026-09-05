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

- **状態**: **保留**(2026-09-05 発見。`autonomy-policy.md` §2「Skeleton の statement の修理」に該当)
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

- **状態**: **保留**(2026-09-06 発見。`autonomy-policy.md` §2 に該当)
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

- **状態**: **保留**(2026-09-06 発見)
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
- **決定**: —
