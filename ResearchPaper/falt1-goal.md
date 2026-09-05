# [Falt1] トラックのゴール(2026-09-04 設定・確定、Chapter I は 2026-09-04 完了)

対象: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299
(物理46ページ、`papers.json`に短縮タグ`Falt1`で登記済み、`0_Source`にPDF/txtあり)。

## 0. ゴール(現在地)

> **Falt1 Chapter I §1 2/2, §2 4/4, §3 2/2, §4 5/5 —— Chapter I 13/13 達成(2026-09-04)**
>
> 残り: Chapter II §1 0/2, §2 0/4, §3 0/2,
> Chapter III §1 0/3, §2 0/2, §3 0/3, §4 0/1 —— 合計 0/17。
> 全体では 13/30。

★★**正直な注記(過大申告を避けるため)**: Chapter I 13項目はすべて
`Interface`(posit、証明していない)+`Skeleton`(.src・.nonvacuous 付き)
+ 260dpi目視確認+ HTML逐語アンカーまで揃い、`node tools/check.mjs`・
`node tools/graph.mjs --owner Falt1`(条なし13・項目13・sorry0)を通る。
ただし:
- §4.4(Theorem 4.4)は[LocProP]形式化時に作った**(i)(iv)のみの縮小形**
  のまま(原文は(i)-(vii)まである)。`.needs`に理由を明記済み。
- §4.2・§4.5は原文が持つ添字(i・M・R1,R2 の分解等)に関する**全称量化を
  代表的な1インスタンスに簡略化**して posit した(`Interface/Falt1/
  GaloisCohomologyCompute.lean` 冒頭に記録)。
- 「almost mathematics」の核心である almost étale covering・m-torsion
  等は mathlib に対応物が無いため全項目が posit ベース(`sorry` は無いが
  `Found`(証明)ではない)。数学的に**証明した**わけではない。

★★★**2026-09-04、`/goal` で「Chapter I を Found(証明)にする」という
より厳しい条件が設定された——その進捗と実現可能性の評価を追記する。**

## 0.0 現在地サマリ(2026-09-05 セッション終了時点、次回の入口)

`/goal` 条件(§1 2/2・§2 4/4・§3 2/2・§4 5/5、計13項目)に対し、
**Found 2/13**(`Lemma 1.1`・`Definition 2.1`)。以下、13項目それぞれの
状態と、次に何をすべきかを1箇所にまとめる(詳細な経緯は本文の★印
エントリを参照、日付は全て2026-09)。

| 節 | 項目 | 状態 | 次の一手 |
|---|---|---|---|
| §1 | Lemma 1.1 | ✅Found | (完成) |
| §1 | Theorem 1.2 | ブロック | `Wₙ⊗_{Vₙ}Vₙ₊₁`の非正規性(`differentIdeal Wₙ Wₙ₊₁`)の評価が必要。`Vₙ`塔側の技術(Eisenstein多項式・differentIdeal計算)は完成済みだが`Wₙ`側(任意のalmost étale `W`)には転用できない(`differentIdeal_tower_diamond`のdocstring参照)。退化した`W`(非分岐)を使うと`δₙ≡0`は自明に出るが、原文の「任意の`W`」という全称量化を証明したことにはならない。 |
| §2 | Definition 2.1 | ✅Found | (完成、`p`が単元でなくても任意のétale/finite/free拡大について成立する一般形まで) |
| §2 | Theorem 2.2 | **一意性✅完成・存在は`p`可除な底待ち** | ★2026-09-05: 入力(remark 2.1(v)、honest 版`hochschild_ext_eq_zero`・almost 版`hochschild_ext_almost_zero`)は完成。**一意性は完全に証明した**(`AlmostLifting.thm_2_2_uniqueness_of_isAlmostEtale`——差が導分になり、`Ω[B⁄A]`のalmost消滅と`C`の`p`捩れ無しから`0`)。**存在**も部品はすべて揃った:第1段`almost_lift_of_isAlmostEtale`(almost射影性からの持ち上げ)・障害のコサイクル性(`obstruction_mem_ker`・`obstruction_identity`)・コバウンダリ表示(`hochschild_H2_almost_coboundary`、明示形)・**Faltings の「`ε`を倍にする」段**(`doubling_multiplicative`)・整合性(`rescale_multiplicative`)。残るのは`ε`族の極限——★これには**`p`可除な底(`m = ∪p^{1/p^k}`)を持つ almost mathematics の層**が要る(下記「構造上の発見」参照)。 |
| §2 | Theorem 2.3 | 同上(同じ変形理論) | `Theorem 2.2`と同じ層を待つ。入力側(`hochschild_ext_almost_zero`は次数を問わず一般)は共通。 |
| §2 | Theorem 2.4(i) | **余核側は完成・核側は1本だけ残る** | ★2026-09-05: `Found/Falt1/AlmostDifferentials.lean`。**余核側は完全に閉じた**——Jacobi–Zariski 完全列(`KaehlerDifferential.range_mapBaseChange`)で余核はちょうど`Ω[B⁄A]`になり、`p^n`が`Ω[B⁄A]`を零化すること(`kaehler_almost_zero`)を条件(iii)の witness から直接示した(`v := w - p^n`と置くと`v∈I`かつ`x·w=0`なので`p^n·x = -(v·x) ∈ I²`——古典的な分離冪等元の1行議論の almost 版)。**核側**は`Algebra.H1Cotangent.exact_δ_mapBaseChange`により「`p^n`が`H1Cotangent A B`を零化する」(= `B`が`A`上 almost formally smooth)にちょうど帰着することまで証明済み(`thm_2_4_i_kernel_of_h1Cotangent`)——残るのはその1本で、これは`Theorem 2.2`(nilpotent ideal に沿った lifting)と同じ深さ。 |
| §2 | Theorem 2.4(ii) | **✅証明済み**(`thm_2_4_ii`、逸脱2件を記録の上) | ★2026-09-05: 4つの部品を全て証明し、原典の仮定の形のまま繋いだ。(1) remark(iii)の trace 恒等式(`remark_iii_trace_identity`)・ノルム適用(`trace_ideal_pow_mem_traceIdeal`)——この2つは`IsAlmostEtaleCovering`仮定の形へ一般化済みで、**`B`自体が étale である必要は無い**(Faltings の設定そのもの)。(2) 群コホモロジー側の transfer を**全次数`i>0`**で(`transfer_groupCohomology_smul_eq_zero`)——以前「mathlibに一般のtransfer定理が無い」と壁として報告していたが、必要なのはrestriction-corestrictionではなく「`Σ_g g(b)=c`なら`c`が`H^i`を零化する」という平均化だけで、coinduced加群への almost split + Shapiro(`groupCohomology.coindIso`)で閉じた。(3) 後半の`M^G/tr_G(M)`(`transfer_invariants_mem_trace`)。(4) **可換環の Galois trace 公式**(`trace_eq_sum_of_chr`——CHRの同型`B⊗_AB≅Map(G,B)`から。mathlibの`trace_eq_sum_automorphisms`は体専用なので自作、`mathlib-gap.json`の`falt1-galois-trace-rings`を3/3で解消)と局所化からの descent(`trace_formula_of_localized`)。**残る逸脱2件**: `Module.Free A B`(原典は`B[1/p]`が projective——mathlibの`Algebra.trace`が`Module.Free`要求、既記録)と`IsDedekindDomain`+`Module.IsTorsionFree`(`Ideal.relNorm`使用のため)。非空虚性の対照(実際にGalois被覆になる具体例)は未作成。 |
| §3 | Theorem 3.1 | ブロック | `Theorem 2.2`-`2.4`の結果を直接使う(§2の壁がそのまま継承) |
| §3 | Theorem 3.2 | ブロック | 同上 |
| §4 | Theorem 4.1-4.3, 4.5 | ブロック | Galois cohomology・スペクトル系列のalmost退化(§2の壁の上にさらに層がある、§2完成が前提) |

**このセッションで新規に確立した再利用可能な事実**(すべて
`lean/ABC3/Found/Falt1/AlmostEtale.lean`、`lake build`・`node tools/
check.mjs --brief`で検証済み、mathlibに無かった一般定理):
- `elem_unique_of_props`:`Algebra.FormallyUnramified.elem`の一意性。
- `diagonalCompare_elem_eq`:`elem`の局所化に関する自然性(`p`が単元
  でなくても成立)。
- `isAlmostEtaleCovering_of_etale_general`:`Definition 2.1`の
  non-vacuous witness、任意の`p`・任意のétale/finite/free拡大。
- `Tr1map_elem_eq_one`:trace form 非退化性(基底の添字計算で証明)。
- `remark_iii_trace_identity`・`trace_ideal_pow_mem_traceIdeal`:
  `Theorem 2.4(ii)`の証明の計算部分。
- `hochschild_ext_eq_zero`:**remark 2.1(v)honest版**——formally
  unramified拡大のHochschild cohomology(`Ext_{S⊗RS}(S,-)`として
  定式化)が正の次数で消える。以前「Faltings未証明」と誤認していた
  事実を実際に証明で訂正した。
- `diagonalCompare_injective`:`B`が`A`上free・`algebraMap B Bp`が単射
  という2条件だけから`diagonalCompare p`が単射(`Module.Flat`の
  rTensor/lTensor保存性2回 + `IsLocalization.moduleTensorEquiv`)。
- `almost_swap_annihilate`・`almost_swap_augment`・
  `hochSectionOfWitness`・`hochSectionAlmost_augment`:
  `IsAlmostEtaleCovering`のみ(honest な`Algebra.FormallyUnramified
  A B`を要求しない)から、`S`が`S⊗AS`の"almost direct summand"
  (`μ∘s=`「`p^n`との掛け算」)であることの形式化。前回セッションで
  報告していた「swap-annihilationを`B⊗AB`へ降ろせない」という技術的
  障壁は、`diagonalCompare_injective`によって解消済み。
- `ext_smul_eq_zero_of_almost_split`:圏論だけで書ける一般補題——
  `s≫μ=τ•𝟙`(almost split)なら`τ`は`Ext^{k+1}(S,M)`を零化する。
- `hochschild_ext_almost_zero`:**remark 2.1(v)almost版、完成**——
  `IsAlmostEtaleCovering A B p`+`Module.Free A B`+`algebraMap B
  B[1/p]`の単射性から、`p^n`が`Ext^{k+1}_{B⊗_AB}(B,M)`の全ての元を
  零化する。Faltings の remark(v)最終文「for B almost etale over A,
  m annihilates the Hochschild cohomology in positive degrees」
  そのもの。真の非単元`p:=5`・`B:=Fin 2 → ℤ`での非空虚性の対照つき。

- `transfer_H1`・`transfer_H2`(コサイクル水準の明示公式)・
  `transfer_H1_smul_eq_zero`・`transfer_H2_smul_eq_zero`(mathlib の
  `groupCohomology.H1`/`H2`に関する形)・
  **`transfer_groupCohomology_smul_eq_zero`(全次数`i>0`)**・
  `coind_almost_split`・`transfer_invariants_mem_trace`:
  **`Theorem 2.4(ii)`の群コホモロジー的内容**——`Σ_{g∈G}g(b)=c`なら
  `c`が`H^{i+1}(G,M)`を零化する。証明の骨格は remark 2.1(v) almost 版と
  同じ「almost split ⟹ almost 消滅」(今セッションで確立したパターンの
  2度目の適用)。非空虚性の対照として古典的な「`|G|`が正次数の群
  コホモロジーを零化する」が系として出ることを確認済み。
- `trace_pi`・`trace_eq_sum_of_chr`・`trace_formula_of_localized`・
  `thm_2_4_ii`:**可換環の Galois trace 公式**(mathlib の
  `trace_eq_sum_automorphisms`は体専用なので自作、CHR の同型
  `B⊗_AB≅Map(G,B)`から`LinearMap.trace_baseChange`+`trace_pi`で導く)と
  局所化からの descent、そして**`Theorem 2.4(ii)`本体を原典の仮定の形
  (Galois は`B[1/p]`の水準)のまま**証明したもの。
- `kaehler_key`・`kaehler_almost_zero`・`thm_2_4_i_cokernel`・
  `thm_2_4_i_kernel_of_h1Cotangent`(`AlmostDifferentials.lean`):
  **`Theorem 2.4(i)`の余核側**(`p^n`が`Ω[B⁄A]`を零化する)と、
  核側が`H1Cotangent A B`の almost 消滅1本に帰着することの証明。
- `almost_dual_basis`・`almost_projective_factor`・
  `almost_lift_of_surjective`・`almost_lift_of_isAlmostEtale`
  (`AlmostProjective.lean`):remark(iii)が与える almost dual basis から
  「`B`は`A`上 almost 射影的」を出し、**`Theorem 2.2`の第1段**
  (任意の全射に沿って`A`-加群写像が`p^n`倍で持ち上がる)を閉じたもの。
- `hochH2Bil`・`hochH2_identity`・`hochschild_H2_almost_coboundary`・
  `obstruction_mem_ker`・`obstruction_identity`・`rescale_obstruction`・
  `doubling_multiplicative`・`uniqueness_derivation_eq`・
  `rescale_multiplicative`(`HochschildLowDegree.lean`):
  **Hochschild `H²`の almost 消滅をコチェインの明示形で**(remark(v)の
  縮約ホモトピーそのもの)与え、`Theorem 2.2`の障害類の扱い
  (コサイクル性・`ε`を倍にする段・一意性・整合性)を**すべて**閉じたもの。
- `derivation_almost_zero`・`thm_2_2_uniqueness`・
  `thm_2_2_uniqueness_of_isAlmostEtale`(`AlmostLifting.lean`):
  **`Theorem 2.2`の一意性を完全に証明**。

### ★構造上の発見(2026-09-05)——なぜ `2.4` は閉じて `2.2`/`2.3` は閉じないか

本プロジェクトの `Definition 2.1` 条件(iii)は指数を **`n : ℕ`** で取る。
原典は `ε>0` を **`p`冪分母の正有理数**で走らせる——`A` が `p^{1/p^k}` を
含む(perfectoid 的な)底を前提とし、`m = ∪_ε p^ε A` は `(p)` より真に
小さい。この違いが項目を二分する:

- **`Theorem 2.4` は almost な結論**(「`m` が零化する」)なので整数指数の
  枠組みでもそのまま意味を持ち、実際に閉じた。
- **`Theorem 2.2`/`2.3` は honest な結論**(「**一意に**持ち上がる」)であり、
  `p^ε φ_ε` の族から `ε→0` の極限として honest な `φ₀` を得る操作が要る。
  整数指数では `ε ≥ 1` までしか下がらず、**`p` で割る操作が原理的に構成
  できない**。一意性側はこの極限を要しないので閉じた(`AlmostLifting.lean`)。

したがって次に着手すべきは「`Theorem 2.2` の続き」ではなく、
**`p`可除な底(`ϖ : ℕ → A`、`ϖ 0 = p`、`(ϖ (k+1))^q = ϖ k`、
`m = ⨆ₖ (ϖ k)`)を持つ層の設計**である。今セッションで揃えた材料は
その層の上でそのまま再利用できる。

**§2 の残りが1点に収束した(2026-09-05)**: `Theorem 2.2`・`2.3`・
`2.4(i)`の核側は、いずれも**「square-zero(冪零)拡大に沿った almost
lifting」という同一の1本**に帰着する:

- `Theorem 2.2` そのものがその主張(`φ:B→C/I` が almost 一意に持ち上がる)。
- `Theorem 2.3` は almost étale covering 自体の持ち上げ(同じ変形理論)。
- `Theorem 2.4(i)`の核側は `p^n·H1Cotangent A B = 0`
  (= `B` が `A` 上 almost formally smooth)で、これは square-zero 拡大に
  沿った almost lifting の言い換えである。

入力(`hochschild_ext_almost_zero`)も、障害類の取り出しも、
`ε`を倍にする段も、今セッションで**すべて明示的に閉じた**
(`HochschildLowDegree.lean`)。**残るのは`ε`族の極限を取る操作1点**で、
それには上記の「`p`可除な底」の層が要る
(`ResearchPaper/mathlib-gap.json`の`falt1-almost-mathematics`)。
次回セッションはこの層の設計が最優先。

## 0.1 `/goal Falt1 Chapter I Found` の進捗(2026-09-04)

**Lemma 1.1**(§1 1/2): ★★★★★**完成**。`Found/Falt1/KaehlerAux.lean`・
`Found/Falt1/Lemma11.lean` として数学的内容(単射性・余核の長さの等式)
を**完全に sorry無しで証明した**(`falt1MapBaseChangeInjective`・
`falt1CokernelLengthEq`・`lemma_1_1_falt1`)。途中で作った再利用可能な
一般補題(mathlib に無かったもの): `subsingleton_H1Cotangent_self`
(`FormallySmooth R S ⟹ Subsingleton(H1Cotangent R S)`)・
`polynomialKaehlerSplit`(`Ω_{V[T]/Z} ≅ (V[T]⊗_VΩ_V)×V[T]` の直和分解)・
`mapBaseChange_injective_transport`(`mapBaseChange` の単射性を代数同型
に沿って輸送)。`differentIdeal_ne_bot` の全条件を `falt1_
differentIdeal_ne_bot` で導出——`differentIdeal V W ≠ ⊥` は仮定では
なく定理になった。★★★(2026-09-04 完成)`Interface.RamificationSetup`
への正式な差し替えも完了した(`falt1RamificationSetup`)——`Module.
length` の `ℕ∞→ℕ` 有限性は当初の見積り(CRT+局所化不変性)より単純に
出た(`Ring.DimensionLEOne`+`Ideal.krullDimLE_zero_quotient_iff_
forall_minimalPrimes_isMaximal`経由)。`RamificationSetup` の全
フィールドを本物のデータで埋め、`lemma_1_1 (falt1RamificationSetup
w hint hadjoin hw)` が型検査を通ることを確認した——Skeleton の posit
(`RamificationSetup.example`)を初めて実データに置き換えた例。
残る記録事項は「絶対微分の基底 `Z`」を一般化(`Z:=ℤ` を正準に選択)
として証明したことのみ(弱化ではない、`.needs` に記録)。

**Theorem 1.2**(§1 2/2): ★★**未着手。260dpi で物理p.5(印字p.258)を
直接読み(pdftotext の OCR は数式記号を激しく壊すため、上記 .txt 行
223-235 は信頼できない——実際に画像で読み直した、2026-09-04)、
実現可能性を評価した結果、Lemma 1.1 より明確に難しいことが分かった。**

- 設定: `V=V_0⊂V_1⊂V_2⊂...` という拡大の列で、`Ω_{V_n/V_{n-1}}`
  (★添字は `V_{n-1}` であって `V` ではない——OCR 行 223-235 の誤読を
  訂正)が `(V_n/pV_n)^{d+1}` を商に持つという条件を満たすもの(典型例:
  `V = W(k)(T_1,...,T_d)` の完備化、`V_n` は `T_i` の `p^n` 乗根と
  `p^{n+1}` 乗根の1のべき根を添加した拡大——perfectoid 系列の原型)。
- 主張: `W_n :=` `V_n⊗_VW` の正規化として、`W_n/V_n` の different
  `δ_n` が `n→∞` で `0` に収束する。
- 証明(原文どおり、正確に読み直した): 塔 `V_n→V_{n+1}→W_{n+1}`
  (と `W_n→W_{n+1}`)から `Ω_{W_n/V_n}⊗_{W_n}W_{n+1} →
  Ω_{W_{n+1}/V_n} → Ω_{W_{n+1}/V_{n+1}} → 0`(完全、これ自体は既存の
  `KaehlerDifferential.exact_mapBaseChange_map` で足りる)。
  **★★新たに要る部品**: 「`Ω_{W_{n+1}/V_n}` が `W_{n+1}/p^αW_{n+1}`
  型の加群 `d+1` 個の直和である」という構造的事実——これは Lemma 1.1 を
  `V_n→V_{n+1}` の**各生成元ごとに**適用して得られるはずだが、原文は
  この構成を「明らか」として省略しており、**単一生成元版の Lemma 1.1
  (今回証明した形)をそのまま複数生成元へ一般化する新規の仕事**になる。
  その後、`Module.length` の完全列に対する加法性・不等式(カーネル・
  コカーネルの長さ評価)を連鎖させ、最終的に `δ_n-δ_{n+1} ≥
  β-(d+1)(δ_n-δ_{n+1})`(`β=min{1,δ_n/(d+1)}`)という漸化不等式から
  `δ_n→0` を導く。

★★★**この構造から、Theorem 1.2 は Lemma 1.1 と同等かそれ以上の
新規インフラ(複数生成元版 Lemma 1.1・`Module.length` の不等式操作の
連鎖・「非常に分岐した」列の定義そのものの形式化)を要求する**——
Lemma 1.1 1件の形式化がこのセッションの大半(数十回のラウンド)を
要したことを踏まえると、Theorem 1.2 は**単独でさらに同程度以上の
作業**になる見込み。★急いで低品質な形式化を試みるより、正確に
読み直した上でこう記録する方が誠実だと判断した。

★★★2026-09-04、**「複数生成元版 Lemma 1.1」の核心(2生成元の場合)が
完成した**(`Found/Falt1/KaehlerAux.lean` の `pushoutKaehlerSplit`)。
`Algebra.IsPushout R C D B`(`B` が `C`・`D` の `R` 上の pushout=
テンソル積)のとき、`Ω_{B/R} ≃ₗ[B] (B⊗_CΩ_{C/R}) × (B⊗_DΩ_{D/R})`
(`mapBaseChange R C B` が単射という条件つき)——`polynomialKaehlerSplit`
と全く同じ「`Function.Exact.splitSurjectiveEquiv` + 明示的セクション」
のパターンで構成できた。**鍵になった発見**: `tensorKaehlerEquivBase`
(片方の脚基点)ではなく **`tensorKaehlerEquiv`**(もう1つの、`B⊗_D
Ω_{D/R} ≃ₗ[B] Ω_{B/C}` という形の mathlib 既製品)の方が正しい道具
だった——`map R C B B ∘ mapBaseChange R D B = tensorKaehlerEquiv
R C D B` という自然性(`tensorKaehlerEquiv_eq_map_mapBaseChange`、
新規に証明)が、`polynomialEquiv` の代わりのセクションを直接与えた。

★残る作業(次のセッションはここから):
1. **`d+1` 個への一般化**(今は2個のみ)——`pushoutKaehlerSplit` を
   `d+1` 回の pushout の反復に適用する帰納法(構造自体は難しくないはず)。
2. ★★★★**2026-09-04、`mapBaseChange R C B` の単射性の条件——
   完成した**(`mapBaseChange_injective_adjoinRoot_tensor`)。下の
   「3つ目(AdjoinRoot基底変換)」の経路が閉じた: `tensorPolynomialAlgEquiv
   (C⊗[R]Polynomial R ≃ₐ[C] Polynomial C)` → `adjoinRootTensorEquiv
   (C⊗[R]AdjoinRoot g ≃ₐ[C] AdjoinRoot(g.map(algebraMap R C)))` →
   `mapBaseChange_injective_of_nzd` + `mapBaseChange_injective_transport`
   で結ぶ。**`B = C⊗[R]AdjoinRoot g` の場合、`pushoutKaehlerSplit` の
   `hinj` は Lemma 1.1 とまったく同じ形の非零因子条件 1 つ
   (`g.map(algebraMap R C)` の微分が `AdjoinRoot(g.map(algebraMap R C))`
   で非零因子)に帰着する。** lake build 済み・sorry 皆無確認済み
   (commit `f45fa30b`)。残る技術的な下ごしらえ: `PUnit` の宇宙を
   `Type`(宇宙0)に固定する必要があった(`R C` を `Type*` にした瞬間
   宇宙メタ変数 `?u` が解決不能になる実測——`tools/lean-idioms.md` へ
   追記の価値あり)。
   ★★2026-09-04(この項目は上で解消済み、経緯として残す)経路を検討・
   部分的に実測した: `Algebra.
   Extension.cotangentComplex_injective_iff P [FormallySmooth R P.Ring]
   : Function.Injective P.cotangentComplex ↔ Subsingleton(H1Cotangent
   R A)`(mathlib)を `P := Extension.ofSurjective(algebraMap
   (Polynomial V)(AdjoinRoot g))` に適用すれば `[FormallySmooth V
   (Polynomial V)]`(既にある)から `Subsingleton(H1Cotangent V
   (AdjoinRoot g)) ↔ Function.Injective P.cotangentComplex` が出る
   はず——問題は **`P.cotangentComplex` が既存の `kerCotangentToTensor`
   と型レベルで(`rfl` でも)一致しないこと**(`P.Cotangent` は定義上
   `P.ker.Cotangent` で `P.ker` も定義上 `RingHom.ker(algebraMap
   P.Ring S)` のはずだが、`Extension.ofSurjective` の内部実装が
   `algebraMap` を直接使わないため、期待通りに簡約されない)。両者を
   **関数として**(`⇑` レベルで)結びつける橋渡しが未解決——`Extension`
   の内部実装(`RingTheory/Extension/Cotangent/Basic.lean` 等)を
   もう少し詳しく調べる必要がある。

   ★★★2026-09-04、**別の(より具体的な)経路も検討した**: `B = C⊗_R
   D`(`D=AdjoinRoot g`)は `C⊗_RPoly R ≅ Poly C`(`MvPolynomial.
   algebraTensorAlgEquiv`、単変数の場合に相当)+ テンソル積が商と
   可換という事実を経由すれば **`B ≅ AdjoinRoot(g.map(algebraMap R C))`
   ——つまり `B` は `C` 上の monogenic 拡大そのもの**になるはず。
   これが成り立てば `mapBaseChange R C B` の単射性は
   `mapBaseChange_injective_of_nzd`(`V:=C`・`f:=g.map(algebraMap R C)`)
   を**そのまま再利用できる**(H1Cotangent 経由の抽象論より遥かに
   具体的)——ただし `C⊗_RPoly R ≅ Poly C` の同型と「テンソル積は
   商と可換」を貼り合わせて `B ≅ AdjoinRoot(...)` を確立する作業
   自体は未着手(`MvPolynomial.algebraTensorAlgEquiv` を単変数
   `Polynomial` へ specialize する下ごしらえが要る)。**3つの経路
   (H1Cotangent+平坦底変換・Extension.cotangentComplex 橋渡し・
   AdjoinRoot 基底変換の同一視)を検討し、いずれも即座には閉じ
   なかった**——この時点ではここで打ち切った。
   →★★★★続きは同日内、上の「完成した」を見よ:3つ目
   (AdjoinRoot 基底変換)を実際に最後まで実行し、閉じた。
★★★★★2026-09-04、**Theorem 1.2 の証明全文(物理p.5=印字p.258)を
260dpi で precisely 読み直した**(OCR ではなく目視)。以下、逐語に近い
形で構造を記録する(今後の形式化がここから始められるように)。

> Proof. Replacing V by some V_n we may assume that all W_n are integral
> domains. Look at the sequence of maps Ω_{W_n/V_n}⊗_{W_n}W_{n+1} →
> Ω_{W_{n+1}/V_n} → Ω_{W_{n+1}/V_{n+1}}. The kernel of the second map
> contains Ω_{V_{n+1}/V_n}⊗_{V_{n+1}}W_{n+1}, which has (W_{n+1}/pW_{n+1})^{d+1}
> as quotient. As Ω_{V_{n+1}/V_n}⊗_{V_{n+1}}W_{n+1} is the direct sum of
> d+1 modules of the form W_{n+1}/p^αW_{n+1}, the kernel of the second
> map contains the kernel of multiplication by p on Ω_{W_n/V_n}⊗_{W_n}W_{n+1},
> and hence the composition of the two maps annihilates the kernel by
> p-multiplication on Ω_{W_n/V_n}⊗_{W_n}W_{n+1}. If we denote the different
> of W_n over V_n by δ_n, this kernel has length at least equal to that
> of W_{n+1}/p^βW_{n+1}, where β = min{1, δ_n/(d+1)}. Also it is clear
> that p^{δ_n-δ_{n+1}} annihilates W_{n+1}/(W_n⊗_{V_n}V_{n+1}), and hence
> the cokernel of the composition of the two maps is annihilated by
> p^{δ_n-δ_{n+1}}. So its length is at most equal to that of W_{n+1}
> divided by the (d+1)st power of this. We derive that
> δ_n-δ_{n+1} ≥ β-(d+1)(δ_n-δ_{n+1}). So if δ_n ≥ d+1, then
> δ_{n+1} ≤ δ_n-1/(d+2), and otherwise δ_{n+1} ≤ (1-1/((d+1)(d+2)))δ_n.
> In any case δ_n → 0 for n→∞. □

直後の「典型例」段落: `V = W(k)(T_1,...,T_d)` の完備化(`W(k)` は完全体
`k` の Witt ベクトル)、`V_n` は「`T_i` の `p^n` 乗根 + `1` の `p^{n+1}`
乗根」で生成される拡大——つまり **`V_{n+1}/V_n` は `d+1` 個の
単生成拡大(各 `T_i` の追加の `p` 乗根 1 個 + 冪根 1 個)の
**同時**添加**(逐次ではない)。これは `pushoutKaehlerSplit` の
`d+1` 項版が直接与える「`Ω_{V_{n+1}/V_n}` は `d+1` 個の巡回加群の
直和」という一文と**厳密に一致する**——★上の「3つの経路」で
確立した `mapBaseChange_injective_adjoinRoot_tensor` と合わせて、
`pushoutKaehlerSplit` は当初の想定通り Theorem 1.2 の証明の**まさに
その箇所**(第3文)で使われるべき道具であることが、憶測ではなく
原文の記述と照合して確認できた。

★以上を踏まえ、残りの形式化作業を(既存の(3)(4)を置き換えて)
以下のように精緻化する:

3a. **`pushoutKaehlerSplit` の `d+1` 項化**——`Fin (d+1)` 添字族
    `D : Fin (d+1) → Type*`(各 `R`-代数、`Algebra.IsPushout` の
    反復合成で `B := D 0 ⊗_R D 1 ⊗_R ⋯ ⊗_R D d` を構成)に対し
    `Ω_{B/R} ≃ₗ[B] ⊕ᵢ (B ⊗_{D i} Ω_{D i/R})`。★この段階では純粋に
    Kähler微分の代数——`V_n`族そのものへの依存は無い、汎用インフラ
    として先に作れる。

    ★★★★2026-09-04、**`d=2`(3項)の場合を最後まで実行し、帰納の
    1ステップが機械的に閉じることを実証した**(`pushoutKaehlerSplit3`、
    commit `d8341aab`)。`pushoutKaehlerSplit` を2回・
    `TensorProduct.AlgebraTensorModule.cancelBaseChange`(「反復底変換
    =直接底変換」)・`TensorProduct.prodRight`(テンソル積と直積の
    可換性)の組み合わせで構成できた——`d+1` 項への一般化もこの2つの
    道具の反復のみで閉じる見込みが実証された。`Algebra C0 B`・
    `Algebra C1 B` とその `IsScalarTower` は(`pushoutKaehlerSplit`
    自身と同じ流儀で)**仮定として明示的に渡す**必要がある——`B1`・`B`
    を具体的な `TensorProduct` として書いた特殊化版は、宣言の**型
    そのもの**が `Algebra.TensorProduct.rightAlgebra`(`abbrev`、
    tools/lean-idioms.md #23)を要求するため別途書けないことも実測した
    (使用箇所ごとに `letI := Algebra.TensorProduct.rightAlgebra` を
    先に置いてから直接呼ぶ必要がある)。
    次のセッションへ: `Fin (d+1)` 一般の帰納は、**型族を再帰的に
    構成する**部分(`B : ℕ → Type*` を `TensorProduct` で再帰定義し、
    各段の `CommRing`/`Algebra` インスタンスを同時に構成する)が
    Lean4 では地雷("recursive Type-valued def with dependent
    instances")——`structure RAlg (R) where carrier; [ring];[alg]`
    のような Σ 束ねを検討すること(未着手)。

    ★★2026-09-04、この一般化への着手を**意図的に保留した**——下の 3b
    の難度を精密に再評価した結果、`pushoutKaehlerSplit` の `d+1` 項化
    (純粋な Lean エンジニアリング)を仕上げても、3b(独立の length
    計算、Lemma 1.1 同等の数学的深さ)が閉じない限り Theorem 1.2 全体
    には接続しない、と判断したため。3a を一般化する労力の投資対効果は
    3b の見通しが立ってから再評価する。
3b. **`Ω_{V_{n+1}/V_n}⊗W_{n+1}` が `d+1` 個の巡回加群の直和という
    事実そのものを `V_n` 族の定義に接続**(3a の適用)。

    ★★★★★2026-09-04、**この項目の難度を精密に再評価した**(原文
    証明の中段落を再度 260dpi で読み直し、独立に再構成を試みた結果)。
    `pushoutKaehlerSplit3` が与えるのは「直和分解の**存在**」だが、
    `delta_tendsto_zero` が要求する `hrec`(漸化不等式)を実際に導く
    には、原文の中段落——「第1の写像(`Ω_{Wₙ/Vₙ}⊗Wₙ₊₁→Ωₙ₊₁/Vₙ`、
    タワー `Vₙ→Wₙ→Wₙ₊₁` の base change)と第2の写像(`Ωₙ₊₁/Vₙ→
    Ωₙ₊₁/Vₙ₊₁`、`KaehlerDifferential.exact_mapBaseChange_map` により
    その kernel はタワー `Vₙ→Vₙ₊₁→Wₙ₊₁` の base change の**像**と一致)
    の**合成**の kernel・cokernel を、`δₙ`(`Wₙ/Vₙ` の different の
    長さ)を使って評価する」という**独立の length 計算**が要る。
    特に「`p^{δₙ-δₙ₊₁}` が `Wₙ₊₁/(Wₙ⊗_{Vₙ}Vₙ₊₁)` を零化する」という
    一文は、`Wₙ⊗_{Vₙ}Vₙ₊₁` と(その正規化である)`Wₙ₊₁` の差を
    differentで測る、半局所環の conductor 理論に相当する事実で、
    **Lemma 1.1 自体と同等かそれ以上の深さの独立した数学的主張**
    であり、`pushoutKaehlerSplit3` を「適用するだけ」では閉じない
    ——単なる工学的な貼り合わせ作業ではなく、**新しい数学的補題
    (半局所 Dedekind 的環の塔での conductor/different の劣加法性)
    を要する**、という結論に達した。★安易な誤形式化を避けるため、
    この評価が固まるまで 3b への Lean 実装は保留する。

    ★★★★★★★2026-09-04、**上記の等式を実際に Lean で構築・証明した**
    (`differentIdeal_tower_diamond`、`Found/Falt1/KaehlerAux.lean`、
    commit `937de066`、sorry無し)——`Algebra.IsSeparable` の引数を
    `FractionRing.liftAlgebra` に明示的に固定しないと2回目の適用で
    instance 探索が食い違うことが判明し(tools/lean-idioms.md #23と
    同種)、それを回避して完成させた。下記の分析(未知数
    `differentIdeal Wₙ Wₙ₊₁` の評価が核心の困難として残る)は変わらず
    有効——これは Theorem 1.2 3b の**代数的骨格**として独立した
    検証済みの一歩である。

    ★★★★★★2026-09-04(発見時のメモ、上で実装完了): mathlib に
    **`differentIdeal` の塔の公式**が既にある
    (`differentIdeal_eq_differentIdeal_mul_differentIdeal`、
    `RingTheory/DedekindDomain/Different.lean`)——`A→B→C` の塔に対し
    `differentIdeal A C = differentIdeal B C * (differentIdeal A B).
    map(algebraMap B C)`(`A` 整閉整域・`B,C` Dedekind・有限・捩れ無し・
    `Algebra.IsSeparable (FractionRing A)(FractionRing C)` が仮定)。
    これを `Vₙ→Vₙ₊₁→Wₙ₊₁` と `Vₙ→Wₙ→Wₙ₊₁` の2通りの塔に適用すると
    `δₙ₊₁ · (differentIdeal Vₙ Vₙ₊₁).map(...) = (differentIdeal Wₙ
    Wₙ₊₁) · δₙ.map(...)` という**厳密な等式**(原文の不等式より強い形)
    が直接得られる——ただし残る未知数 `differentIdeal Wₙ Wₙ₊₁`
    (`Wₙ⊗_{Vₙ}Vₙ₊₁` が正規化前の状態からどれだけ `Wₙ₊₁` と異なるかを
    測る量)を評価する必要は結局残り、これが原文の「kernel・cokernel
    の length 評価」の**まさにその部分**に相当する——つまりこの公式は
    3b の**別の切り口**を与えるが、核心の困難(非正規性の評価)は
    避けられない。ただし「厳密な等式」から出発できる分、不等式の
    往復より見通しが良くなる可能性がある。次のセッションでの
    検討候補として記録する(まだ Lean での検証・使用はしていない)。
★★★★★★★★2026-09-04、**決定的な発見**: `ResearchPaper/0_Source/Brinon
Conrad - CMI Summer School Notes on p-adic Hodge Theory.txt` の
Exercise 13.7.4 が、Faltings のこの証明を**教科書レベルの精密さで
step-by-step に書き下している**(Brinon-Conrad はスタンフォード/
プリンストンで定番の CMI サマースクールノート)。原文(Faltings 自身
の1988年の論文)の圧縮された記述より遥かに読みやすく、以下の疑問点
が解消した:

- (1) `Ω¹_{B1/A}` が `d+1` 個の元で生成される、という事実は
  **「second fundamental exact sequence + Nakayama」で出す**(具体的な
  手法が明示されている——`pushoutKaehlerSplit` 系ではなく、もっと
  直接的な「単項生成の完全列」の議論)。
- (2) **`length_B(Ω¹_{B/A}) = length_B(B/p^{δ_{B/A}}B)`**(`δ_{B/A}`
  は different `D_{B/A} = p^{δ_{B/A}}B` の指数)——これは
  `Found/Falt1/Lemma11.lean` の `falt1CokernelLengthEq` と**本質的に
  同じ計算**(Lemma 1.1 そのもの)であることが確認できた。
- (4) `ker(b) ⊇ ker(p 倍写像)` の正確な意味と根拠が判明:
  「(1)(=d+1個生成)+ elementary divisor theorem」から出る——
  `Ω¹_{An+1/An}⊗Bn+1` が(有限生成 torsion 加群の構造定理により)
  `d+1` 個の巡回加群の直和に分解でき、各巡回加群の位数が `p^{何か}`
  なので、"合成 `b∘a` の kernel は `p` 倍写像の kernel を含む"
  という主張の正確な機構が判明した(★私が「よく分からない」と記録
  していた箇所)。
- (5) **鍵となる未解決の箇所の出典が判明**: 「using the definition of
  the discriminant, show `p^{δn-δn+1}·(Bn⊗_{An}An+1) ⊆ Bn+1`」——
  これは**判別式(discriminant)の塔での乗法性 + 非極大な order の
  conductor が判別式の比較から評価できる、という古典的整数論**
  (conductor-discriminant の関係)を使う、と明記されている。
  `differentIdeal_tower_diamond`(different の塔の公式)ではなく
  **discriminant の塔の公式**が本来の道具だった可能性が高い——
  次のセッションでは mathlib の `Algebra.discr`/`NumberField.discr`
  周りの資産を調査すること(未着手)。

★結論(更新): 核心の困難(非正規性の評価)は依然として**独立の
古典的整数論(conductor-discriminant の関係)を要する**という
これまでの評価は変わらないが、**その具体的な形(discriminant の
塔の乗法性)が判明した**ことで、次にどの mathlib 資産
(`Algebra.discr` 系)を調査すべきかが明確になった。これは
「思ったより深い」から「深いが具体的に何を調べればよいかが分かった」
への前進である。

★★★★★★★★★2026-09-04、続けて mathlib を調査し、**まさにこの
「conductor via discriminant/微分」の道具そのもの**を発見した:
`conductor_mul_differentIdeal`(`RingTheory/DedekindDomain/
Different.lean`)—— `(conductor A x) * differentIdeal A B =
Ideal.span {aeval x (derivative (minpoly A x))}`(`x : B` が
`Frac(B)/Frac(A)` を生成する場合)。これは Lemma 1.1 で使った
`differentIdeal_eq_span_derivative`(`conductor = (1)` の特殊形、
`B` が `A` 上 monogenic な場合)の**一般化**であり、`conductor A x`
が `A[x]`(= 私の `T_{n+1} = Wₙ⊗_{Vₙ}Vₙ₊₁` に相当、`x` が「典型例」
の生成元(`p`乗根等)なら monogenic)が `B`(= `Wₙ₊₁`)からどれだけ
離れているかを**まさに測る量**。つまり:

`differentIdeal Wₙ Wₙ₊₁ = Ideal.span{aeval x (deriv(minpoly Wₙ x))}
/ (conductor Wₙ x)`(概念的に——実際には `conductor * different =
span{...}` の形で使う)

という式で、目標の `Jₙ := differentIdeal Wₙ Wₙ₊₁` が
`conductor Wₙ x` 経由で書ける。★「典型例」(`V_{n+1}/V_n` が
`p`乗根の同時添加)なら `x` の具体的な選び方も明確。

★★★★★★★★★★2026-09-04、**実際に組み合わせて完成させた**
(`cancel_conductor_delta`、`Found/Falt1/KaehlerAux.lean`、commit
`bf0a3ea7`、sorry無し)。`differentIdeal_tower_diamond` と
`conductor_mul_differentIdeal` の2つの等式を掛け合わせるだけで `Jₙ`
が消え、**`conductor(Wₙ,x) · δₙ₊₁ = δₙ.map(...)` というきれいな
関係式**が(「`conductor Wₙ x` と `differentIdeal Vₙ Vₙ₊₁` の
base change が一致する」という `hspan_eq` 仮定のもとで)出ることを
確認した。Dedekind整域の「0でないイデアルは消去可能」という性質
だけで閉じる、純粋に代数的な議論。

★残る作業: (a) `hspan_eq`(3c=「典型例」で `x` の最小多項式が base
change と両立することの確認)、(b) `conductor(Wₙ,x)` 自体の下限
評価(`p^{何か}` を含むことの証明、Brinon-Conrad step (5) 後半に
対応)、(c) この関係式(掛け算の形)から `delta_tendsto_zero` が
要求する `hrec`(引き算の不等式の形)への変換。

★★★★★★★★★★★2026-09-04、(c) に実際に着手し、**「(a)(b)より軽い」
という見立てが誤りだったと判明した**。`length(R/(I·J)) =
length(R/I) + length(R/J)`(Dedekind整域、`I,J` 0でない、必ずしも
coprimeでない一般の場合)という長さの加法性が要るが、これは
mathlib に無い(coprimeの場合の `Ideal.quotientMulEquivQuotientProd`
はあるが一般の場合は無い)。証明の筋は分かる——完全列
`I/(IJ) → R/(IJ) → R/I → 0`(`Module.length_eq_add_of_exact` が
使える)+ `I/(IJ) ≅ R/J`(可逆イデアル `I` によるテンソルが
length を変えないという事実、局所化して「可逆イデアルは各素点で
自明」に帰着させる議論)——だが後半の同型は**局所化と length の
両立性**という、これも mathlib に無さそうな一般論を要する。
★このセッションで繰り返し遭遇したパターン(「見た目は軽いはずが
掘ると独立の一般補題に当たる」)がここでも再現した——安易に
「軽い」と決めつけず、正直に記録する。次のセッションへ:
`length(R/(I·J)) = length(R/I)+length(R/J)` を独立の補題として
先に確立すること(Dedekind整域論のごく基本的な事実のはずだが、
mathlib での正確な組み立て方は未確認)。

    ★★2026-09-04、この長さの加法性を埋める手がかりを追加調査した。
    `Ideal.relNorm : Ideal S →*₀ Ideal R`(`RingTheory/Ideal/Norm/
    RelNorm.lean`、**任意の Dedekind 底環 `R`** に対する相対ノルム、
    `→*₀` として bundled——`ℤ` 限定の `Ideal.absNorm` と違い一般の塔
    に使える)は型として乗法的性なので `relNorm(IJ)=relNorm(I)*
    relNorm(J)` はタダで手に入る。ただしこれは `S` 上の長さではなく
    `R` の**イデアルへ**写す不変量なので、`length_S(S/I)` を直接
    与えるには別途「`length_S(S/I) = length_R(R/relNorm I)`」的な
    比較公式が要り、それ自体が `R` 側で**同じ加法性**を要求する
    (循環)——この経路は今回は閉じなかった。次のセッションへ:
    「単一の素イデアル `P` に対し `length_S(S/P^n) = n·length_S(S/P)`」
    を最初の一段として直接攻める(`TensorProduct.tensorQuotEquivQuotSMul`
    で `P/P^{n+1} ≅ S/P` を示し帰納する筋、`P` が可逆という事実を
    経由)方が、`relNorm`/`absNorm` 経由より見通しが良いかもしれない。

    ★★★★★★★★★★★★2026-09-04、**核心の道具を発見した**:
    `RingTheory/PicardGroup.lean` の `Module.Invertible`(可逆加群、
    Picard群論のための typeclass、`contractLeft` の全単射性で定義)に、
    まさに必要な「テンソルは完全性を保つ」補題群がある:
    `Module.Invertible.rTensor_injective_iff`・`_surjective_iff`・
    `_bijective_iff`(`[Module.Invertible R M]` のとき
    `f.rTensor M` の単射性・全射性・全単射性は `f` 自身のそれと同値)。
    これは「不変(=可逆)加群でテンソルすると完全列が完全列のまま」
    という、`Module.length_eq_add_of_exact` と組み合わせれば
    `length(M⊗N)=length(N)` を**帰納法で**導出できるはずの、
    まさに探していた一般論——mathlib に無いと思っていたが、
    「Picard群」という一見無関係な名前のファイルに実在した
    (`measure-mathlib-before-skeleton` の教訓が今回も再現: 概念名
    ではなく数学的内容で探すべきだった)。

    ★残るブリッジ: `IsDedekindDomain S`(の nonzero イデアル `I`)から
    `Module.Invertible S ↥I` インスタンスを得る経路——
    `isDedekindDomain_iff_isDedekindDomainInv`(`RingTheory/
    DedekindDomain/Ideal/Basic.lean`)経由で `I * I⁻¹ = 1`
    (`FractionalIdeal` として)は出るが、これを `Module.Invertible`
    の `contractLeft` 全単射性の形へ変換する配線はまだ未着手。

    ★★★★★★★★★★★★★2026-09-04、**この一般化の代わりに、単項イデアルの
    場合に限定して完成させた**(`length_quotient_span_singleton_mul`、
    `Found/Falt1/KaehlerAux.lean`、commit `25221786`、sorry無し)。
    `Module.Invertible`(可逆加群・局所化)の一般論は一切不要で、
    既存の道具(`Module.length_eq_add_of_exact`・`Submodule.
    quotientQuotientEquivQuotient`・`LinearMap.quotKerEquivRange`)
    の組み合わせだけで、`a` が非零因子のとき
    `length(R/((a)·J)) = length(R/(a)) + length(R/J)` を証明できた。
    多くの実際の場面(`differentIdeal_eq_span_derivative`・
    `conductor_mul_differentIdeal` 等、単一の生成元の微分で書ける
    場合)ではイデアルが単項なので、3b(c) の核心 gap をこの特殊形で
    埋められる見込みが高い。★これで3b(c)は実質的に解決した——
    残るのは `cancel_conductor_delta` の各イデアル(`conductor Wₙ x`
    等)が実際に単項であることの確認(3a と同様、「典型例」の具体形
    次第)。

    ★★★★★★★★★★★★★★2026-09-04、残る (a)(b) の依存関係を精査した。
    (a)(`minpoly Wₙ x = (minpoly Vₙ の生成元).map(...)`という
    base change 整合性)は `minpoly.eq_of_irreducible_of_monic`
    (mathlib、体上でのみ成立)を経由するのが自然だが、これには
    「`f.map(algebraMap Vₙ Wₙ)` が `Frac(Wₙ)` 上既約のまま」という
    非自明な仮定が要り、これは**抽象的に証明できる一般命題ではなく、
    具体的な `Vₙ`・`Wₙ` の構成(3c)に依存する事実**である(一般には
    偽になりうる——`Wₙ` が `Vₙ` の拡大の情報を「知りすぎている」と
    既約性が壊れる)。(b)(conductor の下限評価)も同様に、具体的な
    ramification の値(`p` の分岐指数等)を要する。
    ★★結論: **(a)(b) は 3b の中の独立した補題としてではなく、3c
    (「非常に分岐した」`V_n` 族の具体的構成)の一部として初めて
    意味を持つ**——3b(c) のように「3c 抜きで先に済ませられる」
    部分ではなかった。次に着手すべきは 3c そのもの
    (3a の `d+1` 項化・3b(a)(b) いずれも、3c が具体的な対象を
    与えて初めて検証可能になる)。

    ★★★★★★★★★★★★★★★★★2026-09-04、`cancel_conductor_delta` を
    `delta_tendsto_zero` の `hrec` へ実際に接続しようとして、
    もう1つの繋ぎ目(「次数によるスケーリング」: `Wₙ₊₁` が `Wₙ` 上
    有限自由のとき `length_{Wₙ₊₁}(Wₙ₊₁/I.map(...)) = [Wₙ₊₁:Wₙ]·
    length_{Wₙ}(Wₙ/I)` という形の等式)が要ると気づいたが、
    **これは誤りだと検算で判明した**——記録しておく価値がある失敗:
    `Wₙ₊₁/Wₙ` が(残余体拡大のみの)**不分岐**次数 `n` の場合を検算
    すると、`Wₙ₊₁/I·Wₙ₊₁` はそれ自体が(残余体の)体になり
    `Wₙ₊₁`-加群としての長さは **1**(不分岐がどれだけ高次でも)。
    一方 `length_{Wₙ}(Wₙ/I)` も体なので **1**。つまり `[Wₙ₊₁:Wₙ]=n`
    倍にはならない——**正しい係数は次数 `n` 全体ではなく分岐指数
    `e` のみのはず**(残余次数分は両辺の「1加群あたりの大きさ」に
    暗黙に吸収されるため)。★誤った式のまま次に進まず、ここで
    立ち止まって検算したことを記録する——このセッションを通じて
    「見た目の一般化を急がず、具体例(今回は不分岐の場合)で検算する」
    姿勢が繰り返し役に立っている。

    ★★★★★★★★★★★★★★★★★★2026-09-04、**正しい形を完成させた**
    (`length_map_pow_of_ramificationIdx`、`Found/Falt1/KaehlerAux.lean`、
    commit `f3e5c53d`、sorry無し)。`Wₙ₊₁` が DVR で、その極大イデアルが
    `Wₙ` のイデアル `P` の像として `𝔪^e`(`e`=分岐指数)と書けるとき、
    `P^k` の像で割った長さは `e·k`——`IsDiscreteValuationRing.length_
    quotient_pow_maximalIdeal` から直接従う、短い証明だった。これで
    「次数スケーリング」の繋ぎ目も(局所=DVRの場合に)閉じた。

3c. **「非常に分岐した」`V_n` の族**そのものの形式化(具体例: `p^n`
    乗根と `1` のべき根を添加する塔、上の「典型例」段落)——まだ
    手つかず。抽象的な族の公理化(`Ω_{V_n/V_{n-1}}` が
    `(V_n/pV_n)^{d+1}` を商に持つ、という条件だけを構造として
    持たせる)なら 3a・3b より軽い可能性がある——ただし 3b が独立の
    length 計算を要すると分かった今、3c だけを先に進めても
    Theorem 1.2 全体の完成には直結しない(3b を経由しない限り
    `delta_tendsto_zero` の `hrec` に接続できない)。

    ★★★★★★★★★★★★★★★2026-09-04、3b(c) 完成の勢いで 3c 着手を検討した。
    フル(Witt ベクトル)の「典型例」ではなく、**最小の非自明な事例
    (`d=0`、単一の一様化元の逐次 `p` 乗根: `V_{n+1} := V_n[X]/(X^p-π_n)`
    、`π_n` は `V_n` の一様化元)**を試みた。この場合
    `X^p-π_n` は `V_n` 上 **Eisenstein**(`Polynomial.IsEisensteinAt`、
    mathlib にあり、`.irreducible` も既製品)で、これが 3b(a) の
    「既約性が保たれる」を**自動的に**満たす、という良い性質がある。
    ★ただし `AdjoinRoot(X^p-π_n)` が実際に整閉(`V_{n+1}` の役を
    果たせる、classical に「Eisenstein ⟹ 極大整環」)であることの
    直接の mathlib 定理はまだ見つけていない——`conductor_mul_
    differentIdeal` 経由で「Eisenstein なら conductor=(1)」を
    示す形で回り道できる可能性はある(未検証)。時間の制約で
    ここで打ち切ったが、**d=0 の最小事例は3a(生成元1個、既存の
    Lemma 1.1 インフラをそのまま使える)・3b(a)(Eisenstein の既約性
    で自動)の両方を回避できる有望な入口**と判断する——次のセッション
    はここから始めるのが良い。

    ★★★★★★★★★★★★★★★★2026-09-04、**「Eisenstein ⟹ 整閉」の道具を
    さらに調査した**: `RingTheory/Polynomial/Eisenstein/IsIntegral.lean`
    に、まさにこの classical な証明の**部品**(`mem_adjoin_of_smul_
    prime_pow_smul_of_minpoly_isEisensteinAt`:「`p^n•z` が
    `adjoin R{gen}` に入れば `z` も入る」)が既にあることを発見した
    ——これは教科書の証明(Eisenstein 多項式の根で生成される拡大が
    極大整環であることの証明)の**核心の帰納法ステップ**そのもの。
    ただし**「任意の整な `z` に対しある `n` で `p^n•z∈adjoin`」という
    近似の議論(教科書証明のもう半分)へつなぐ具体的な組み立ては
    まだ行っていない**——これも classical だが、mathlib に「Eisenstein
    ⟹ 整閉」という完成形の定理としては見つからなかった(`totally
    ramified`という語もmathlibに無い)。★これで3cの技術的な下ごしらえ
    (Eisenstein 既約性・整閉性それぞれの部品の所在)はかなり明確に
    なった——次のセッションは `mem_adjoin_of_smul_prime_pow_smul_of_
    minpoly_isEisensteinAt` を軸に「整閉性」の証明を完成させることから
    始めるとよい。

    ★★★★★★★★★★★★★★★★★★★2026-09-04、**より見通しの良い別証明の筋
    を発見した**(未着手・未検証): `mem_adjoin_of_smul_prime_pow_smul_
    of_minpoly_isEisensteinAt` 経由の一般論(length・局所化理論を要し
    このセッションで繰り返し壁に当たった領域)ではなく、**付値による
    直接計算**の方が自己完結している可能性が高い。`A` を(一様化元
    `π` を持つ)DVR、`f:=X^n-π`(Eisenstein)、`θ:=` `AdjoinRoot f` の
    根とすると: (i) `Frac(B)` 上へ `A` の付値 `v` を延長すると
    `v(θ)=1/n·v(π)` になる(`f` が Eisenstein なので強制される)、
    (ii) `y=c_0+c_1θ+...+c_{n-1}θ^{n-1}`(`c_i∈Frac(A)`)の付値は
    `min_i(v(c_i)+i/n)`(**`i/n` が `i=0,...,n-1` で mod 1 相異なる
    ので相殺が起きない**、超距離不等式の非退化版)、(iii) これから
    直接「`y` が整 ⟺ 全ての `c_i∈A`」が従う——つまり
    `AdjoinRoot(X^n-π) = B` の整閉性そのものが、**局所化やlengthの
    一般論を経由せず、付値の直接計算だけで**出る。★mathlib に
    「有限拡大への付値の延長」・「相異なる付値を持つ項の和は最小値を
    与える」という部品が使える形であるかは未確認——次のセッションの
    最有力候補としてこちらを先に調べる価値がある(`mem_adjoin_of_...`
    経由より見通しが良い可能性)。

    ★★2026-09-04、上の「付値の延長」に使える可能性のある mathlib の
    部品の**在り処だけ**確認した(中身は未検証):`Valuation.
    HasExtension`(`RingTheory/Valuation/Extension.lean`)・
    `Valuation.extendToLocalization`(`RingTheory/Valuation/
    ExtendToLocalization.lean`)——`KaehlerAux.lean` の現在の import
    には含まれず、別途 import してから中身を調べる必要がある。

    ★★★2026-09-04、中身を確認し、**`extendToLocalization` はこの
    用途に使えないと判明した**——名前通り「局所化」への延長専用
    (`S ≤ v.supp.primeCompl`・`IsLocalization S B` が要る)であり、
    `B := AdjoinRoot(f)`(`A` の局所化ではなく、次数 `n` の有限拡大)
    には適用できない。`Valuation.HasExtension` は延長の**存在**では
    なく**両立性**(`vR` が `vA` の comap と同値)を述べる Prop で、
    構成自体は別途必要。★「有限次拡大への付値の延長」という、まさに
    必要な一般論の所在はまだ見つかっていない——次のセッションでは
    `IsDedekindDomain.HeightOneSpectrum.valuation` 経由(`B` 自身が
    Dedekind であることを**先に**示してから、その height-one spectrum
    の付値として `θ` の付値を得る、という順序の入れ替え)を検討する
    余地があるが、これは「`B` が Dedekind」という結論を仮定して使う
    ことになり得るため、循環を避ける組み立てが必要。

    ★★★★★★★★★★★★★★★★★★★★2026-09-04、**戦略そのものを転換する
    重要な発見をした**: そもそも `AdjoinRoot(f)` が「既に」整閉である
    ことを Eisenstein 性から証明する必要は**無かった**。mathlib には
    `integralClosure.isDedekindDomain_fractionRing`(`RingTheory/
    DedekindDomain/IntegralClosure.lean`、**instance**)があり、
    `A` が Dedekind・`L` が `FractionRing A` 上有限分離拡大なら
    `integralClosure A L` は**自動的に** Dedekind になる
    (`FractionRing A` を明示的に使うことが鍵——一般の `K` with
    `IsFractionRing A K` では instance が発火しない、という
    tools/lean-idioms.md #23 と同種の罠を実測で確認・回避した)。
    **`Vₙ₊₁ := integralClosure Vₙ L` と定義すれば、Eisenstein 性・
    付値の延長・整閉性の証明が一切不要になる**——`w`(`f` の根)が
    `Vₙ₊₁` を「生成し切っているか」は問わなくてよい。なぜなら
    `conductor_mul_differentIdeal`(**既に完成済み**)が、その
    conductor を未知数のまま追跡できるよう作られているから——
    `cancel_conductor_delta` はまさにこの「`w` が生成し切らない
    かもしれない」場合を扱うために作った道具だった。★これで
    3c の技術的な核心の壁(整閉性の証明)を**丸ごと迂回できる**
    ことが判明した。次のセッションは `integralClosure.isDedekindDomain_
    fractionRing` を使って実際に `V_1 := integralClosure V_0 L`
    (`L:=FractionRing(AdjoinRoot(X^p-π))`)を組み立て、Falt1
    RamificationSetup 系のインフラ(`falt1_isFractionRing` 等)を
    再利用して具体的な最初の1段を完成させることから始めるとよい。

    ★★★★★★★★★★★★★★★★★★★★★2026-09-04、**上記の戦略転換を実装し、
    「非常に分岐した」`V_n` 族の1段の構成が sorry 無く完成した**
    (commit `9998b221`)。`isDedekindDomain_integralClosure_adjoinRoot_
    X_pow_sub_C`: `V0` を Dedekind、`π` を(平方でない)素元とし
    `f := X^n - C π`(`V0[X]`上)とすると、`integralClosure V0
    (AdjoinRoot (f.map (algebraMap V0 (FractionRing V0))))` が
    Dedekind であることを証明。組み立ての鍵:
    - `f` の `V0[X]` 上の Eisenstein 既約性(既存の
      `eisenstein_X_pow_sub_C_irreducible_map`)を**Gauss の補題**
      (`Polynomial.Monic.irreducible_iff_irreducible_map_fraction_map`
      、`IsIntegrallyClosed R` だけで足り **UFD は不要**——Dedekind
      整域は自動的に整閉なのでこれで賄える)で `FractionRing(V0)[X]`
      上の既約性に転送する。これで「`AdjoinRoot(f)` が整域か」という
      (UFD でない一般の Dedekind 整域上では自明でない)問いを完全に
      回避できた。
    - 分離性は `Polynomial.separable_X_pow_sub_C` を直接使う(自作の
      Bezout 係数計算は不要と判明し破棄)。
    - `AdjoinRoot(f_K)`(`f_K:=f.map(...)`、`FractionRing(V0)` 上)が
      有限次分離拡大であることを、`Module.Finite`
      (`Polynomial.Monic.finite_adjoinRoot`)と新設の補題
      `algIsSeparable_adjoinRoot_of_separable`(単既約 monic 分離
      多項式の `AdjoinRoot` は分離拡大——
      `IntermediateField.isSeparable_adjoin_simple_iff_isSeparable`
      + `adjoin_root_eq_top` + `topEquiv` を貼り合わせて構成)から得る。
    - 最後に `integralClosure.isDedekindDomain_fractionRing`
      (instance)を `infer_instance` で発火させるだけで結論が出る。

    `#print axioms` で両定理とも `[propext, Classical.choice,
    Quot.sound]` のみ(sorry 無し)を確認済み、`lake build`
    2656/2656 成功、`node tools/check.mjs --brief` は新規リグレッション
    無し(既存 13 件の NG のみ)。★これで item 3c の「`V_n` 族の各段が
    Dedekind であること」という核心のボトルネックが、**`n=0→1` の
    1段について**完全に解消した。残る作業: (a) この構成を `V_0`
    自身にも再帰的に適用して `V_n`(`n` 段目)の族を作る(`induction`
    で回すだけのはずだが、各段で「一様化元」を選び直す必要があり
    未着手)、(b) `cancel_conductor_delta` の `hspan_eq` 仮定
    (「conductor が単項で書ける」)をこの具体的な `V_1` 構成に対して
    実際に検証する(3b との接続、まだ着手していない)、(c) `Ω¹_{V1/V0}`
    が `pushoutKaehlerSplit` 系の枠組みに乗ることの確認(生成元の個数
    `d+1` との対応)。

    ★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**単項生成性
    (monogenicity)も完成した**(commit `3e1f9798`)——`falt1CokernelLengthEq`
    が要求する `hadjoin : Algebra.adjoin V ({w}:Set W) = ⊤` を、上の
    `V_1 := integralClosure V0 L` に対して実際に示した。決定的な発見:
    **同じリポジトリの別の論文トラック**(`ABC3.Found.GenEll.
    TameRamification`、Mochizuki [GenEll] の馴分岐の形式化)に、
    まさに必要な部品が既に完成していた——
    `mem_adjoin_of_pow_smul_of_isEisensteinAt`(Eisenstein なら
    `π^k•z∈adjoin` から `z∈adjoin`、mathlib の `mem_adjoin_of_smul_
    prime_smul_of_minpoly_isEisensteinAt` を `k` 回反復する帰納)と
    `exists_smul_mem_adjoin_powerBasis`(`L` の任意の元はある `0≠d∈R`
    倍すれば `adjoin` に入る、`IsLocalization.exist_integer_multiples`
    経由)。**DVR** 上では任意の非零 `d` が `π^k×単元`
    (`IsDiscreteValuationRing.eq_unit_mul_pow_irreducible`)と分解
    できるので、この2つを橋渡しするだけで「任意の整な元は `adjoin`
    に入る」(`adjoin_eq_top_of_isEisensteinAt`)、つまり `adjoin R
    {gen} = integralClosure R L` が出た。★これは教科書の「Eisenstein
    拡大は単項生成(totally ramified ⟹ monogenic)」という古典定理
    そのものだが、`AdjoinRoot(f)` の整閉性を直接証明する迂回路を
    通らずに到達できた——`Subalgebra.map_injective` +
    `AlgHom.map_adjoin_singleton` で「`L` の中での等式」を「`W` 自身の
    中での `⊤` 生成」に変換し(`adjoin_eq_top_in_integralClosure_of_
    isEisensteinAt`)、`minpoly V0 w = f`(`X^n-π` そのもの、`minpoly.
    isIntegrallyClosed_eq_field_fractions'` + `Polynomial.map_injective`
    で `K` 上の minpoly `f_K` から `V0` 上の minpoly へ降ろす)を経由
    して具体例へ組み立てた
    (`adjoin_eq_top_integralClosure_adjoinRoot_X_pow_sub_C`)。
    `#print axioms` で全定理 sorry 無しを確認、`lake build`
    2837/2837 成功、`node tools/check.mjs --brief` は新規リグレッション
    無し。★これで **`V_1` の1段について、Dedekind であること
    (前回)と単項生成であること(今回)の両方が揃った**——
    `falt1CokernelLengthEq` を実際に適用する準備が整った。次回は
    上記の残る作業 (a)(b)(c) のうち、特に (b)(`hspan_eq` の検証、
    `differentIdeal V0 V1 = span{f'(w)}` を `differentIdeal_eq_span_
    derivative`(Lemma 1.1 自身の道具)で計算し、`cancel_conductor_
    delta` の等式と突き合わせる)に取り組むのが、3b と 3c を実際に
    接続する最短経路と見込む。

    ★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**(b) を実行し、
    `differentIdeal V0 V1` の具体式が完成した**(`differentIdeal_eq_
    span_of_adjoinRoot_X_pow_sub_C`)。`differentIdeal_eq_span_
    derivative`(Lemma 1.1 自身の道具、`w` が単項生成なら `differentIdeal
    V W = span{f'(w)}`)と、前段の monogenicity を貼り合わせるだけで、
    `differentIdeal V0 V1 = span{n·w^{n-1}}`(`w`=根の像、`n·X^{n-1}`
    は `X^n-π` の微分)という**古典的な Eisenstein 拡大の判別式の公式
    そのもの**が閉じた。★配管上の教訓(新パターン、tools/lean-idioms.md
    に追記の価値あり): `differentIdeal A B` は `[IsDedekindDomain B]`・
    `[Module.IsTorsionFree A B]` を引数に持つため、これらを証明内部で
    `infer_instance` で揃えようとすると、**結論(goal)自体が同じ2つの
    instance を要求する**ことで先に走る instance 探索の失敗が
    キャッシュされ、後から `haveI` で正しく揃えても失敗し続けるという
    罠に嵌った(`set` 起因ではない、新しい失敗形)——解決策はこれら2つを
    **定理の引数として明示的に要求する**ことで goal 自体の instance
    探索を迂回すること。

    ★★委員会メモ(作業の帰属について): この commit は**共有
    worktree で並行して動いていた別セッション**(Lubin-Tate/pGC の
    冪級数評価に取り組んでいたセッション)が `git commit -a`(または
    同等の非pathspec制限のコミット)を実行した際に、**このセッションが
    ステージ済みだった `KaehlerAux.lean` の変更を巻き込んで一緒に
    コミットした**——commit hash は `88879620`(コミットメッセージは
    「冪級数を Λ_n の元で実際に評価する」という無関係な内容だが、
    diff には本節で述べた `differentIdeal_eq_span_of_adjoinRoot_
    X_pow_sub_C`(112行)が含まれている、実測で確認済み)。sorry 無し・
    `lake build` 成功・push 済みであることは確認済みなので**内容自体は
    無事**——今後 `git log -- lean/ABC3/Found/Falt1/KaehlerAux.lean`
    でこの完成を辿る際は `88879620` を見ること。

    ★これで item 3c の V_1 について、**Dedekind・単項生成・
    differentIdeal の具体値**の3点セットが揃った。残る作業:
    `cancel_conductor_delta` の `hspan_eq` 仮定(`spanDeriv = Idiff`)を
    この具体式(`Idiff = Ideal.map (algebraMap V1 W1) (differentIdeal
    V0 V1) = Ideal.map(...) (span{n·w^{n-1}})`)に対して実際に
    代入・検証する接続作業(3b と 3c を本当に繋ぐ最後の1ピース)、
    および上記(a)(V_n 再帰族)・(c)(`pushoutKaehlerSplit` との接続)。

    ★★★★★★2026-09-04、**`hspan_eq` 接続の具体的な設計図が見えた**
    (未着手・次回の最有力候補)。`differentIdeal_tower_diamond` の
    `Wₙ→Wₙ₊₁` 側を、`Vₙ→Vₙ₊₁`(=今回の `V0→V1`)と**同じ多項式
    `X^n-π` を `Wₙ` へ base change したもの**で作れば
    (`Wₙ₊₁ := integralClosure Wₙ (AdjoinRoot((X^n-C π).map(V0→Wₙ→
    Frac(Wₙ))))`)、`hspan_eq` に必要な「`x`(`Wₙ₊₁` の生成元)の
    `algebraMap V1 Wₙ₊₁` による `w` の像との一致」は、**今セッション
    の前半で作った `adjoinRootTensorEquiv`**(`C⊗[R]AdjoinRoot g
    ≃ₐ[C] AdjoinRoot(g.map(algebraMap R C))`、`pushoutKaehlerSplit`
    のために構築したもの)を `R:=V0`・`C:=Frac(Wₙ)`・`g:=fK` で適用
    すれば**ほぼ自動的に従う**——`w⊗1` の像が `AdjoinRoot` 側の
    `root(fK.map(...))` に対応するという、`adjoinRootTensorEquiv` の
    定義性質そのもの。つまり **`Wₙ₊₁` 側の全ての情報
    (Dedekind・単項生成・differentIdeal の値)も、`differentIdeal_
    eq_span_of_adjoinRoot_X_pow_sub_C` 系の3定理を**そのまま
    `V0` の代わりに `Wₙ` に適用するだけ**で再利用できる**(3定理は
    すでに一般の `[IsDiscreteValuationRing V0]` に対して書いてある
    ので、`Wₙ` が DVR でありさえすれば即座に適用できる)。
    残るのは(i)`Wₙ` が DVR で `algebraMap V0 Wₙ π` が引き続き
    prime・非自乗であること(`Wₙ/V0` が今考えている素での不分岐、
    という自然な仮定)を明示的な仮説として立てること、
    (ii)`adjoinRootTensorEquiv` の像が `integralClosure` の元として
    振る舞うことの橋渡し(`IsIntegral` が base change で保たれる、
    という一般に易しいはずの事実)、の2点——(ii)は次回最初に
    着手するのに適した、比較的小さい補題だと見込む。

    ★★★★★★★★★2026-09-04、**base change 写像そのものが完成した**
    (`algHomAdjoinRootOfCompat`・`algHomAdjoinRootOfCompat_root`、
    commit `7f09476c`、他セッションに巻き込まれてコミット——上と
    同じ状況、内容は sorry 無し・push 済み確認済み)。当初
    `adjoinRootTensorEquiv`(テンソル積経由)を使う計画だったが、
    実装時により直接的な経路を発見した: `AdjoinRoot.map`(mathlib、
    係数の ring hom `φ0:K→+*K'` から `AdjoinRoot fK→+*AdjoinRoot
    (fK.map φ0)` を作る)を `AlgHom.mk'` で `V0`-線形性
    (`φ0` が `V0` 上の algebraMap と両立するという1つの仮定
    `hφ0` から、`IsScalarTower.algebraMap_apply` +
    `AdjoinRoot.algebraMap_eq` + `AdjoinRoot.map_of` で従う)を
    確認するだけで足りた——`adjoinRootTensorEquiv` も
    `AdjoinRoot.mapAlgHom`(mathlib、同じ底環上のみ)も不要だった。
    `algHomAdjoinRootOfCompat_root` で「根は根に写る」ことも確認
    済み(`w` の像が `Wₙ₊₁` 側の `x` と一致することの鍵)。

    ★残る作業(上記(i)相当、未着手): `φ0 := IsFractionRing.map
    (algebraMap V0 Wₙ の単射性)`(mathlib、`FractionRing V0 →+*
    FractionRing Wₙ`)を実際に構成し、`fK.map φ0` が「`Wₙ` を底環に
    直接構成した `X^n-π'`(`π':=algebraMap V0 Wₙ π`)の Eisenstein
    多項式」と一致すること(`π'` が引き続き prime・非自乗という
    仮定のもとで)を示し、`differentIdeal_tower_diamond`・
    `conductor_mul_differentIdeal` のフル仮説束(`IsDedekindDomain
    Wₙ`・`Module.Finite Wₙ Wₙ₊₁`・`IsScalarTower V0 Wₙ Wₙ₊₁` 等)を
    実際に組み立てて `cancel_conductor_delta` を適用する——これが
    3b・3c を実際に閉じる最終組み立てとして残っている。

    ★★★2026-09-04、**`φ0` 自体も完成した**(`fractionRingMapOfInjective`・
    `fractionRingMapOfInjective_algebraMap`、commit `fb2396c6`)——
    `IsFractionRing.map`(mathlib)+ `IsLocalization.map_eq` +
    `IsScalarTower.algebraMap_apply` で、`algebraMap V0 Wₙ` が単射
    という1つの仮定だけから `φ0:FractionRing V0→+*FractionRing Wₙ` と
    その `V0` 上の algebraMap との両立を得た。これで
    `algHomAdjoinRootOfCompat`(前段)に実際に渡す `φ0`・`hφ0` の
    両方が揃った——**base change 写像の構成は完全に完成**、残る
    作業は上記の「`π'` が引き続き Eisenstein であることの仮定の
    もとでの多項式の一致」と「フル仮説束の組み立て」のみ。

    ★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**驚きの単純化と、
    3c→3b接続の核心が完成した**(commit `2ec2215c`・`865622c8`・
    `f43ae887`)。実装を進める中で、`adjoinRootTensorEquiv`(テンソル積
    経由)や `algHomAdjoinRootOfCompat`(フラクション体経由)より
    **遥かに単純な経路**を発見した:

    1. **`falt1AdjoinRootEquivIntegralClosure`**: `V_1 := integralClosure
       V0 (AdjoinRoot(X^n-π の base change))` は、実は**base change
       する前の** `AdjoinRoot(X^n-π)`(`V0[X]` 上そのままの多項式の商)
       と `V0`-代数として**同型**——`minpoly.equivAdjoin` + monogenicity
       (`adjoin_eq_integralClosure_of_isEisensteinAt`)を貼り合わせる
       だけ。★配線上の罠(記録): `PowerBasis` の暗黙引数を伴う定理を
       裸の `w` に対して呼ぶと高階単一化が `whnf` タイムアウトする
       ——常に `PB.gen` の形を経由すること。
    2. **`algHomAdjoinRootOfCompat'`**: (1)のおかげで、`Vₙ→Vₙ₊₁` の
       橋渡しは **`V0[X]→Wₙ[X]` の底環 base change だけ**(`AdjoinRoot.
       map` + `AlgHom.mk'`)で済むと判明——フラクション体・
       `IsFractionRing.map` は一切不要だった。
    3. **`falt1BaseChangeAlgHom`**: 上記2つを**3段合成**
       (`V1≃ₐ[V0]AdjoinRoot f → AdjoinRoot f→ₐ[V0]AdjoinRoot g
       (`g:=f.map(algebraMap V0 Wₙ)`、多項式として `X^n-Cπ'` に等しい)
       → AdjoinRoot g≃ₐ[Wₙ]Wₙ₊₁`、最後の等式を `AlgHom.restrictScalars
       V0` で `V0`-線形とみなして繋ぐ)して、**`V1 →ₐ[V0] Wₙ₊₁` の
       実際の `AlgHom` を構成した**——`differentIdeal_tower_diamond`
       が要求する `Algebra Vn1 Wn1` インスタンスがこれで手に入る。

    `#print axioms` で全定理 sorry 無し確認、`lake build` 2837/2837
    成功、`node tools/check.mjs --brief` は新規リグレッション無し。
    ★★これで **3c→3b 接続の最大の技術的な壁(「同じ生成元を共有する
    塔をどう構成するか」)を実際に越えた**。残る作業:
    (a) `falt1BaseChangeAlgHom` の下で生成元の対応(`ψ(w)=x` 相当)を
    明示的に確認すること、(b) `Module.Finite Wₙ Wₙ₊₁`・
    `IsScalarTower V0 Wₙ Wₙ₊₁`・`IsScalarTower V0 V1 Wₙ₊₁` 等の残りの
    instance を(これも `falt1AdjoinRootEquivIntegralClosure` 経由で
    「`Wₙ₊₁` も結局 `AdjoinRoot` そのもの」という単純化が使えるはず)
    整備すること、(c) `differentIdeal_tower_diamond` の `hsep`
    (`Algebra.IsSeparable (FractionRing V0)(FractionRing Wₙ₊₁)`)を
    示すこと、(d) 最終的に `conductor_mul_differentIdeal` +
    `cancel_conductor_delta` を実際に適用すること——これで item 3c は
    技術的な核心を越え、残りは「組み立て」の段階に入った。

    ★★★2026-09-04、**(a) に着手し、Lean の配線上の教訓を得た**
    (未完成——`falt1AdjoinRootEquivIntegralClosure` 自体は既に完成・
    push 済みで、以下は「その値を `AdjoinRoot.root f` で評価した式を
    後から取り出す」という**別の**作業)。`unfold` + `simp only
    [eq_mpr_eq_cast, eq_mp_eq_cast, cast_eq]` で `rw […] at e1` が
    生んだ**入れ子の cast の大半は消せる**(外側の2つ、`set f`・
    `set fK` に由来するもの)——ただし `hEqmin`(`minpoly V0 PB.gen =
    f`)由来の**内側の cast**(`AlgEquiv` という関数型そのものへの
    cast)だけは `cast_eq`・`AlgEquiv.trans_apply`・`Subtype.ext_iff`
    のどれでも綺麗に剥がせず、`subst`/`generalize` も(`f` が `set`
    由来の let で「motive is not type correct」)効かなかった——
    各試行が **50〜75秒**かかる中で複数回失敗し、時間対効果が悪化した
    ため一旦中断。

    ★★推奨する迂回路(次回、これを最初に試すこと): **後から `e1
    (AdjoinRoot.root f) = w` を証明しようとするのではなく、`w` を
    最初から `w := e1 (AdjoinRoot.root f)` と定義すればこの等式は
    `rfl` になる**——`differentIdeal_eq_span_of_adjoinRoot_X_pow_
    sub_C` 側の `w`(`⟨PB.gen, hBint⟩` として個別に構成したもの)を、
    `falt1AdjoinRootEquivIntegralClosure` から得られる `w` に**差し
    替える**方向で両者を統一すれば、cast を剥がす作業そのものが
    不要になる可能性が高い。★この種の「opaque な tactic-mode `def`
    を後から `unfold` して cast まみれの中身を分析する」パターンは
    時間対効果が悪いと判明した——tools/lean-idioms.md #29 に記録した。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**この推奨迂回路を
    実行し、生成元の対応が完成した**(`falt1GeneratorPackage`・
    `falt1BaseChangeAlgHom_generator_correspondence`、commit
    `cc34f441`)。予想通り `w := e1 (root f)` と最初から定義する
    ことで、以前 `PowerBasis`・Eisenstein の詳細を経由していた
    「整性・単項生成性・`minpoly=f`」の3点セット
    (`falt1GeneratorPackage`)が、**同型の自然性だけ**
    (`IsIntegral.map`・`AlgHom.map_adjoin_singleton`・
    `minpoly.isIntegrallyClosed_dvd`+`Irreducible.associated_of_dvd`)
    から出ることを確認した——`isDedekindDomain_...`・`adjoin_eq_top_
    ...` 系の元の証明より遥かに短い代替ルート。さらに `falt1BaseChangeAlgHom`
    という名前を経由して `w↦x` を証明しようとする(`unfold` して中身を
    調べる)アプローチは何度試みても50〜75秒ずつタイムアウトし続けた
    ため、**同じ構成をその場で直接書き下す**(named def を経由しない)
    ことで解決した——`algHomAdjoinRootOfCompat'_cast_root`(`g` を
    自由変数のまま保つのが鍵、`set`束縛された変数に対する `subst`
    は「invalid equality proof」で失敗する)を独立補題として切り出し、
    それを直接使うことで `hstep1root` が一瞬(1秒未満)で通った。

    ★★これで **item 3c → 3b 接続に必要な技術的ピースが全て揃った**
    (Dedekind・単項生成・differentIdeal の具体式・base change 写像・
    生成元対応)。残る作業は純粋な「組み立て」段階: (a) `Module.Finite
    Wₙ Wₙ₊₁`・`IsScalarTower V0 Wₙ Wₙ₊₁`・`IsScalarTower V0 V1 Wₙ₊₁`
    等の残りの instance を整備する(これも `falt1AdjoinRootEquivIntegralClosure`
    経由で「`Wₙ₊₁` も結局 `AdjoinRoot` そのもの」という単純化が使える
    はず、`Module.Finite` は `Polynomial.Monic.finite_adjoinRoot` から
    ほぼ自明)、(b) `differentIdeal_tower_diamond` の `hsep`
    (`Algebra.IsSeparable (FractionRing V0)(FractionRing Wₙ₊₁)`)を
    示す、(c) `conductor_mul_differentIdeal` を `Wₙ→Wₙ₊₁` に適用して
    `Jₙ`(`differentIdeal Wₙ Wₙ₊₁`)を conductor で書く、(d) 最終的に
    `cancel_conductor_delta` を適用して `condWnx * deltaN1 =
    deltaNmapped` を得る——ここまで来れば 3b・3c は技術的には閉じ、
    残るのは `delta_tendsto_zero` の `hrec` への変換(3b(c)、既存の
    `length_quotient_span_singleton_mul` 等が使えるはず)と、V_n の
    **再帰族**(item 3c 本体、`V_1` の1段からの帰納)の構築のみ。

    ★★2026-09-04、`differentIdeal_tower_diamond` の instance 整備
    (上記(a))に着手した(commit `45ecea4c`)。**`Module.Finite Wₙ
    Wₙ₊₁` を完成**(`falt1ModuleFiniteWnWn1`)——`AdjoinRoot.powerBasis'`
    (**体でなく一般の可換環上でも** monic なら `PowerBasis` が取れる、
    mathlib の既製品、`Polynomial.Monic.finite_adjoinRoot` の体限定
    版より広い)+`falt1AdjoinRootEquivIntegralClosure` を貼り合わせる
    だけで閉じた。★また **`Algebra V0 Wₙ₊₁`(合成 `V0→Wₙ→FractionRing
    Wₙ→AdjoinRoot gK`)と `IsScalarTower V0 Wₙ Wₙ₊₁` は instance 探索
    だけで自動的に見つかる**ことを実測で確認した(明示的な構成は不要)。
    ★残る instance(未着手): `Module.IsTorsionFree Wₙ Wₙ₊₁`
    (`falt1BaseChangeAlgHom_generator_correspondence` の証明内で
    使った `Module.IsTorsionFree V0 (AdjoinRoot g)` の議論をほぼ
    そのまま流用できる見込み)、`Algebra V1 Wₙ₊₁` を instance として
    登録する(`falt1BaseChangeAlgHom.toRingHom.toAlgebra`)、
    `IsScalarTower V0 V1 Wₙ₊₁`(`falt1BaseChangeAlgHom` が `→ₐ[V0]`で
    あることから `ψ.commutes` 経由でほぼ自動のはず)、そして `Wₙ`
    自身が `V0` 上有限・捩れ無しであること(`Wₙ` がどんな `V0`-拡大
    として与えられるかに依存する、呼び出し側の仮定として要る)。

    ★★2026-09-04、**`Module.IsTorsionFree Wₙ Wₙ₊₁` も完成した**
    (`falt1ModuleIsTorsionFreeWnWn1`、commit `addad6e1`)——予想通り
    `differentIdeal_eq_span_of_adjoinRoot_X_pow_sub_C` の中で使った
    `Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot fK))`
    の議論(`algebraMap Wₙ (AdjoinRoot gK)` の単射性 →
    `IsIntegralClosure.isTorsionFree`)を `Wₙ`・`gK` に対して**そのまま
    繰り返すだけ**で閉じた——新しい発見は無かったが、既存パターンの
    再利用が効くことを確認できた。★これで約20個中、`Module.Finite
    Wₙ Wₙ₊₁`・`Module.IsTorsionFree Wₙ Wₙ₊₁`・`Algebra V0 Wₙ₊₁`・
    `IsScalarTower V0 Wₙ Wₙ₊₁`(自動発見)の4つが揃った。残る
    (未着手・次回): `Algebra V1 Wₙ₊₁` を instance として登録する
    (`falt1BaseChangeAlgHom.toRingHom.toAlgebra`)、`IsScalarTower
    V0 V1 Wₙ₊₁`、`Module.Finite V0 V1`・`Module.IsTorsionFree V0 V1`
    (`V1` 自身、既存の `falt1ModuleFiniteWnWn1`・`falt1ModuleIsTorsionFreeWnWn1`
    と同じパターンを `V0`・`fK` に戻して適用するだけのはず——まだ
    やっていないが機械的)、`Module.Finite V0 Wₙ`・`Module.IsTorsionFree
    V0 Wₙ`(`Wₙ` 自身、呼び出し側の仮定として要る)、`Module.Finite
    V1 Wₙ₊₁`・`Module.IsTorsionFree V1 Wₙ₊₁`(`ψ` 経由の有限性、
    `Module.Finite.trans`的な議論が要りそう)、`Module.Finite V0 Wₙ₊₁`・
    `Module.IsTorsionFree V0 Wₙ₊₁`(`Wₙ→Wₙ₊₁` と `V0→Wₙ` の合成)。

    ★★2026-09-04、**`Module.Finite/IsTorsionFree V0 V1` も完成した**
    (`falt1ModuleFiniteV0V1`・`falt1ModuleIsTorsionFreeV0V1`、commit
    `331b1d9d`)——`Wₙ`側と全く同じ議論を `V0`・`f` に戻すだけ。
    ★★重要な発見: `Module.IsTorsionFree V0 V1` はこれまで**多くの
    定理が instance 前提として明示的に要求していた**が、実は毎回
    この議論で導出可能だったと判明した(既存の証明済み定理群は
    変更せず、影響範囲を数える余裕が無いため追加のみに留めた——
    次回、余裕があれば `[Module.IsTorsionFree V0 (...)]` の前提を
    削って引数を減らすリファクタリングを検討する価値がある)。

    ★★★`Algebra V1 Wₙ₊₁`・`IsScalarTower V0 V1 Wₙ₊₁` も**パターンとしては
    確認できた**(未commit・named theorem としてのパッケージ化は
    断念)——`letI := (falt1BaseChangeAlgHom ...).toRingHom.toAlgebra`
    としてから `IsScalarTower.of_algebraMap_eq` + `ψ.commutes` で
    `IsScalarTower V0 V1 Wₙ₊₁` が出ることを、1つの大きな `example`
    ブロック内で実際に確認した(約19秒、sorry 無し)。★ただしこれを
    **独立した named theorem として `∃`-wrap してパッケージ化しよう
    とすると**(`Algebra` インスタンスを存在型で包んで返す形)、
    `refine ⟨inferInstance, ?_⟩` の instance 探索が
    `«synthesize pending MVars»` でタイムアウトし続けた——
    tools/lean-idioms.md #29 と同種だが、今回は「instance を戻り値と
    して`∃`で包む」という**新しい**失敗パターン。★推奨: この種の
    「これから使う instance を用意する」系の事実は、独立した
    theorem として証明・命名しようとせず、**実際に
    `differentIdeal_tower_diamond` を呼び出す箇所で `letI`/`haveI`
    としてインラインに展開する**のが正攻法——次回はこの方針で
    直接 `differentIdeal_tower_diamond` の呼び出しに挑むとよい
    (単独の instance-provider theorem を目指さない)。

    ★★2026-09-04、**`Module.Finite/IsTorsionFree V0 Wₙ₊₁` も完成した**
    (`falt1ModuleFiniteV0Wn1`・`falt1ModuleIsTorsionFreeV0Wn1`、commit
    `d739bbed`)。`Module.Finite` は `Module.Finite.trans`(`V0→Wₙ→
    Wₙ₊₁` の合成、`Module.Finite V0 Wₙ` を仮定として要求)で機械的。
    `Module.IsTorsionFree` は `algebraMap V0 Wₙ₊₁` の単射性を
    `algebraMap V0 Wₙ`(仮定)・`algebraMap Wₙ Wₙ₊₁`(構成から)の
    合成で示した。★これまで3回繰り返していた「正則元→非零→
    `algebraMap` 単射性→整域での非零元倍は単射」という議論を、
    **汎用補題 `moduleIsTorsionFree_of_injective`**(`T` が整域で
    `algebraMap R T` 単射なら `Module.IsTorsionFree R T`)として
    独立に切り出した——今後同種の場面ではこれを直接呼べばよい。

    ★これで約20個中**8個**(`Module.Finite`/`IsTorsionFree` の
    `(Wₙ,Wₙ₊₁)`・`(V0,V1)`・`(V0,Wₙ₊₁)` の3ペア分、`Algebra V0 Wₙ₊₁`・
    `IsScalarTower V0 Wₙ Wₙ₊₁`)が揃った。残る(未着手):
    `Algebra V1 Wₙ₊₁`・`IsScalarTower V0 V1 Wₙ₊₁`(パターン確認済み、
    `letI`でインライン展開する方針)、`Module.Finite/IsTorsionFree
    V1 Wₙ₊₁`(`ψ:=falt1BaseChangeAlgHom` の単射性を要する、
    まだ証明していない新しい数学的内容)、`hsep`(`Algebra.IsSeparable
    (FractionRing V0)(FractionRing Wₙ₊₁)`)。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**`ψ` の単射性が
    完成した**(`adjoinRoot_map_mk`・`adjoinRoot_map_injective_of_map_eq`・
    `algHomAdjoinRootOfCompat'_injective`・`_cast_injective`・
    `falt1BaseChangeAlgHom_generator_and_injective`、commit `e3ae2da5`)
    ——このセッション後半で唯一残っていた**新しい数学的内容**だった。
    予想(基底 `{1,...,X^{n-1}}` の対応)通りだったが、実際には
    **`PowerBasis` の基底展開を一切経由せず**、mathlib の
    `Polynomial.map_dvd_map`(`x` が monic・係数写像が単射なら
    `x.map f ∣ y.map f ↔ x ∣ y`——多項式の割り算アルゴリズムが monic
    除数に対して任意の可換環上で機能することの帰結)を直接使うだけで
    閉じた——`AdjoinRoot.map` の核が「`p∣r`」という単純な条件に
    帰着することが、`AdjoinRoot.lift_mk`+`Polynomial.eval₂_map`+
    `AdjoinRoot.aeval_eq` で `mk` への作用を計算し、`AdjoinRoot.mk_
    eq_zero` へ橋渡しするだけで示せた。

    ★★続けて **`Module.IsTorsionFree V1 Wₙ₊₁` も即座に完成した**
    (`falt1ModuleIsTorsionFreeV1Wn1`、commit `ce9e4bc2`)——`ψ` の
    単射性を `moduleIsTorsionFree_of_injective`(既存の汎用補題)に
    直接渡すだけ。★これで約20個中**9個**+ `Algebra V1 Wₙ₊₁`・
    `IsScalarTower V0 V1 Wₙ₊₁`(パターン確認済み)が揃った。
    **残るのは `Module.Finite V1 Wₙ₊₁`(単射像の有限生成性の議論、
    `ψ` が単射なのでおそらく `Submodule.fg_of_injective` 的な既製品
    一発で閉じるはず)と `hsep` のみ**——`differentIdeal_tower_diamond`
    の instance 整備は完成が視野に入った。

    ★★★★2026-09-04、**`Module.Finite V1 Wₙ₊₁` も完成した**
    (`falt1ModuleFiniteV1Wn1`、commit `786fcd84`)——予想(単射像の
    有限生成性)とは違う、もっと単純な経路で閉じた:`Module.Finite.
    of_restrictScalars_finite`(`R→A→M` のタワーで `M` が `R` 上有限
    なら `A` 上でも有限、「小さい環上の有限生成は大きい環上でも
    有限生成」という一般に易しい事実)を `R:=V0,A:=V1,M:=Wₙ₊₁` に
    適用するだけ——`IsScalarTower V0 V1 Wₙ₊₁`(`ψ.commutes` から)と
    既存の `Module.Finite V0 Wₙ₊₁` の2つを揃えれば機械的。★`ψ` の
    単射性は**この有限性には不要だった**(`IsTorsionFree` の方でのみ
    必要)——予想が外れたが、より単純な経路が見つかった好例。

    ★★★これで **`differentIdeal_tower_diamond` の instance 整備は
    事実上完成した**——約20個中10個を直接証明し、残り(`Algebra
    Vn1 Wn1`・`Algebra Vn Wn`・`Algebra Vn Vn1`・`IsScalarTower Vn
    Vn1 Wn1`・`IsScalarTower Vn Wn Wn1` 等)は全て確認済みの
    自動発見パターン(`letI`でインライン展開)で揃う。**唯一残る
    のは `hsep`(`Algebra.IsSeparable (FractionRing V0)(FractionRing
    Wₙ₊₁)`)のみ**——ただしこれは技術的に厄介な点が1つ判明した:
    `integralClosure.isFractionRing_of_finite_extension`(mathlib)
    により `IsFractionRing Wₙ₊₁ (AdjoinRoot gK)` は言えるが、これは
    `FractionRing Wₙ₊₁`(**具体的な** `Localization (nonZeroDivisors
    Wₙ₊₁)` という構成)と `AdjoinRoot gK` が**同型**であることしか
    意味せず、**等しい**わけではない——`hsep` を `Algebra.IsSeparable
    (FractionRing V0)(AdjoinRoot gK)`(こちらは `algIsSeparable_
    adjoinRoot_of_separable` 系の既存道具が直接使える形)から
    移送するには、この同型に沿った separability の transport が
    もう1段要る。★また `Wₙ` 自身の `FractionRing V0` 上の分離性
    (`Algebra.IsSeparable (FractionRing V0)(FractionRing Wₙ)`)も
    新たな仮定として要りそう(`Wₙ` が塔のどこから来るかに依存)。
    次回はこの「同型に沿った transport」から始めるとよい——
    tools/lean-idioms.md #29 の教訓(cast を避け、`IsFractionRing`
    の一意性補題 `IsFractionRing.algEquivOfIsFractionRing` 等を
    使って明示的な同型として扱う)が活きる場面のはず。

    ★★★★★★★2026-09-04、**`hsep` の数学的内容そのものは完全に検証
    できた**(未commit・named theorem へのパッケージ化は断念)。
    予想通りの筋(分離性の推移律 + `IsFractionRing` の一意性による
    同型 transport)で完全に閉じることを、**独立した `example` ブロック
    内で sorry 無く実証した**:

    1. `IsLocalization.algEquiv`(mathlib)で `FractionRing Wₙ₊₁ ≃ₐ[Wₙ₊₁]
       AdjoinRoot gK` を得る(`integralClosure.isFractionRing_of_
       finite_extension` が与える `IsFractionRing Wₙ₊₁ (AdjoinRoot gK)`
       から)。
    2. **驚きの発見**: `Algebra (FractionRing V0)(FractionRing Wₙ)`
       が(`fractionRingMapOfInjective` 経由で)一度 `letI` で入って
       さえいれば、`Algebra (FractionRing V0)(AdjoinRoot gK)` と
       `IsScalarTower (FractionRing V0)(FractionRing Wₙ)(AdjoinRoot
       gK)` は**instance 探索だけで自動的に見つかる**(`φ` を手で
       組み立てる必要は無かった)。
    3. `Algebra.IsSeparable.trans`(mathlib、分離性の推移律)で
       `Algebra.IsSeparable (FractionRing V0)(AdjoinRoot gK)` を
       `Algebra.IsSeparable (FractionRing V0)(FractionRing Wₙ)`
       (新たな仮定)+ `Algebra.IsSeparable (FractionRing Wₙ)
       (AdjoinRoot gK)`(`algIsSeparable_adjoinRoot_of_separable`、
       既存)から得る。
    4. `AlgEquiv.ofRingEquiv`(mathlib)で(1)の同型を`FractionRing
       V0` 上の `AlgEquiv` に昇格させ(`e_Wn1.apply_symm_apply` から
       `rfl` 級の可換性証明)、`AlgEquiv.Algebra.isSeparable_iff`
       (mathlib)で(3)を transport すれば `hsep` そのものが出る。

    ★★ただし**この4ステップを1つの named theorem として綺麗に
    パッケージ化しようとすると**(`Algebra.IsSeparable (FractionRing
    V0)(FractionRing Wₙ)` を instance 前提として要求する形にすると)、
    (2)の自動発見が**named theorem の中では効かなくなる**という
    新しい instance-search の脆さに遭遇した——`letI`で導入した
    ローカルな `Algebra` インスタンスと、`hsep` 系の仮定の型が
    参照する(別途 `def` で名指しされた)同じ値が、**instance 探索の
    観点では別物として扱われる**らしく、複数回の書き方の変更でも
    再現し続けた。★これは`falt1BaseChangeAlgHom_generator_and_
    injective`パッケージ化の際に遭遇したものとも異なる**新しい**
    パターン——tools/lean-idioms.md への追記候補(次回、実際に
    追記する前にもう少し原因を特定する価値がある)。

    ★結論: **`hsep` は数学的には解決済み**(4ステップの筋は完全に
    正しく、個別に sorry 無く確認できた)。残るのは Lean の配線上の
    問題(named theorem 化)のみ——次回は、`differentIdeal_tower_
    diamond` を**直接呼び出す箇所**で、このセッションの `example`
    ブロックとまったく同じ順序で `letI`/`haveI` を並べてインライン
    展開する(named theorem を経由しない)方針を試すのが良い
    (`falt1BaseChangeAlgHom_generator_and_injective` 完成時と同じ
    教訓——instance を提供する事実は独立に名前を付けようとせず、
    使う場所で直接組み立てる)。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**実際にこの方針で
    `differentIdeal_tower_diamond` を1つの巨大な `example` で直接
    呼び出してみた——`hsep` 以外の**約19個の instance は全て問題無く
    組み上がる**ことを確認した**(sorry 無し、`hsep` を1つの仮定
    として明示的に受け取る形なら`True`は達成、commit 未——named
    theorem としては未パッケージ)。組み立て順序:
    `isDedekindDomain_integralClosure_adjoinRoot_X_pow_sub_C`(V1・
    Wₙ₊₁ 双方)→ `falt1ModuleIsTorsionFreeV0V1`・`falt1ModuleFiniteV0V1`
    → `falt1ModuleIsTorsionFreeWnWn1`・`falt1ModuleFiniteWnWn1` →
    `falt1ModuleFiniteV0Wn1`・`falt1ModuleIsTorsionFreeV0Wn1` →
    `falt1BaseChangeAlgHom_generator_and_injective` で `ψ` を obtain
    → `letI := ψ.toRingHom.toAlgebra` で `Algebra V1 Wₙ₊₁` →
    `IsScalarTower V0 V1 Wₙ₊₁`(`ψ.commutes`)→
    `moduleIsTorsionFree_of_injective hψinj` →
    `Module.Finite.of_restrictScalars_finite`。**この順序で19個
    全てが素直に通った**——大きな驚きは無く、既存の各ピースが
    そのまま合成できることを確認できただけ、という意味で
    「配線の仕上げ」段階に入ったことの実証。

    ★ただし**`hsep` を最後に接続しようとした瞬間、mathlib 自身の
    `FractionRing.liftAlgebra`(`abbrev`)との diamond
    (tools/lean-idioms.md #23 と全く同じクラスの問題)に遭遇した**:
    `differentIdeal_tower_diamond` が実際に要求する `Algebra
    (FractionRing V0)(FractionRing Wₙ₊₁)` は `FractionRing.
    liftAlgebra V0 (FractionRing Wₙ₊₁)`(`Algebra V0 Wₙ₊₁` から
    自動的に持ち上げる mathlib の既製品)であり、これまでこの
    セッションで組み立ててきた「`AdjoinRoot gK` 経由の手作りの
    instance」とは**別の項**になる。#23 の教訓通り `letI :=
    FractionRing.liftAlgebra V0 (FractionRing Wₙ₊₁)` を明示的に
    呼び、`IsScalarTower V0 (FractionRing V0)(FractionRing Wₙ₊₁)`
    が(このセッションで確認済み、`infer_instance` で自動的に
    見つかる)ことを使って `IsLocalization.ringHom_ext`(mathlib、
    局所化からの2つの準同型が底環上で一致すれば全体で一致する)で
    `hcommutes` を再構成する、という**筋道までは特定した**——
    ただし `algebraMap V0 Wₙ₊₁`(標準、`V0→Wₙ→Wₙ₊₁` 経由)と
    `φ`(`AdjoinRoot gK` へ向かう、`V0→Wₙ→FractionRing Wₙ→AdjoinRoot
    gK` 経由)が同じ元に対して一致する、という**もう1段の
    diamond 解消**(`e_Wn1` を Wₙ の scalar でも可換にする議論)が
    追加で必要と判明し、ここで打ち切った。

    ★結論(正直な評価): `differentIdeal_tower_diamond` の**19/20
    の instance は完全に組み上がることを実証した**——これは
    Theorem 1.2・item 3c にとって非常に大きな前進。**残る `hsep`
    1つ**は、数学的には完全に正しい(既に別の形で証明済み)ものの、
    mathlib の `FractionRing.liftAlgebra` という**具体的な instance
    と一致させる**という、既知のパターン(#23)の**入れ子になった
    版**として立ちはだかっている——次回はこの「二重の diamond
    解消」(`FractionRing.liftAlgebra` との一致 + `V0→Wₙ₊₁` の
    2つの経路の一致)から始めるのが最短距離と見込む。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**`hsep` を完全に
    解決した**(`falt1_hsep_bundled`、`KaehlerAux.lean`、commit
    `cd7a95ec`、`lake build`・`#print axioms` 確認済み・sorry 無し)。
    上で「もう1段の diamond」と呼んでいた懸念——`algebraMap V0 Wₙ₊₁`
    (標準経路)と `φ`(`AdjoinRoot gK` 経由の経路)が一致するか——は
    実は**`rfl` で即座に解消する**ことが判明した(`Wₙ₊₁ =
    integralClosure Wₙ (AdjoinRoot gK)` への代入が `Subalgebra` の
    値そのものなので、`Subtype.val` の展開が両辺で完全に一致する)。
    残る本体は次の二段構成:
    (1) V0-level: `hinjV0Wn1`(`V0→Wₙ₊₁` の単射性、`Polynomial.
    map_dvd_map` 経由の`hinjV0Wn`から)→ `FaithfulSMul V0 (FractionRing
    Wₙ₊₁)`(`faithfulSMul_iff_algebraMap_injective`)→ `letI :=
    FractionRing.liftAlgebra V0 (FractionRing Wₙ₊₁)` で mathlib
    純正の instance を得る(`IsScalarTower V0 (FractionRing
    V0)(FractionRing Wₙ₊₁)` は`infer_instance`で自動的に出る)。
    (2) Wₙ-level でも全く同じパターンを1段繰り返す(`FaithfulSMul Wₙ
    (FractionRing Wₙ₊₁)` → `letI := FractionRing.liftAlgebra Wₙ
    (FractionRing Wₙ₊₁)`)——ここから `FractionRing Wₙ₊₁ ≃ₐ[Wₙ₊₁]
    AdjoinRoot gK`(`IsLocalization.algEquiv`)を`Wₙ`のscalar
    でも可換にする(`e_Wn1.commutes` + 上記`rfl`)ことで
    `Algebra.IsSeparable (FractionRing Wₙ)(FractionRing Wₙ₊₁)` を
    `hsepAdjoin`(`AdjoinRoot gK`側、既存)から移送する
    (`AlgEquiv.ofRingEquiv` + `AlgEquiv.Algebra.isSeparable_iff`)。
    最後に `IsScalarTower (FractionRing V0)(FractionRing
    Wₙ)(FractionRing Wₙ₊₁)`(2つの `FractionRing.liftAlgebra`
    instance 間の両立性、`IsScalarTower.of_algebraMap_eq'` +
    `IsLocalization.ringHom_ext (nonZeroDivisors V0)` で pointwise に
    `IsScalarTower.algebraMap_apply` を繋いで示す)を経て
    `Algebra.IsSeparable.trans` で `(FractionRing V0)→(FractionRing
    Wₙ)→(FractionRing Wₙ₊₁)` を1本に繋いだ。
    ★パッケージング上の教訓(新規): `theorem` の**戻り値の型**
    (ここでは結論)の中で `FractionRing.liftAlgebra V0 (...)` を
    直接使おうとすると、型検査の時点(証明本体に入る前)で
    `FaithfulSMul`・`Fact (Irreducible gK)` 等の instance が要求され、
    シグネチャ側だけでは間に合わない(または `whnf`/`isDefEq` が
    数十秒〜timeoutする)。**`Σ' (inst : Algebra ...), Algebra.
    IsSeparable ... inst` で instance そのものを戻り値として束ね、
    `theorem` ではなく `def`(戻り値が `Type`)にする**ことで、
    シグネチャの型検査が instance 検索を一切要求しない形になり、
    問題が消えた(`lean-idioms.md` に追記候補)。
    これで `differentIdeal_tower_diamond` の**20/20 instance が
    全て揃った**——次回はこの `falt1_hsep_bundled` の `.1`・`.2` を
    上記19-instance の組み立てに接続し、`differentIdeal_tower_diamond`
    を実際に1回呼び出して diamond 等式そのものを得ることが最短距離。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**接続を実行し、
    `differentIdeal_tower_diamond` を Falt1 の `V0→V1`・`Wn→Wn1`
    構成に対して初めて実際に呼び出した**(`falt1_differentIdeal_tower_
    diamond_assembled`、`KaehlerAux.lean`、commit `03a485fa`、
    `lake build`・`#print axioms` 確認済み・sorry 無し)。20個の
    instance(`IsDedekindDomain`×2・`Module.IsTorsionFree`×5・
    `Module.Finite`×5・`Algebra`×5・`IsScalarTower`×2・`hsep`)を
    既存の `falt1Module*` 群 +
    `falt1BaseChangeAlgHom_generator_and_injective`(`ψ` を得る)+
    `falt1_hsep_bundled` から機械的に組み立て、
    `differentIdeal_tower_diamond hsep` を実際に適用して
    ```
    differentIdeal V1 Wn1 * Ideal.map (algebraMap V1 Wn1) (differentIdeal V0 V1)
      = differentIdeal Wn Wn1 * Ideal.map (algebraMap Wn Wn1) (differentIdeal V0 Wn)
    ```
    (`V1 := integralClosure V0 (AdjoinRoot fK)`、`Wn1 := integralClosure
    Wn (AdjoinRoot gK)`、`Algebra V1 Wn1 := ψ.toRingHom.toAlgebra`)を
    型検査済みの項として得た——**これが Theorem 1.2・item 3c の
    中核道具(different の塔の公式)の初めての具体的インスタンス化**。

    ★パッケージング上の限界(記録): この定理の**戻り値**を(呼び出し側が
    再利用しやすいよう)`∃ ψ, letI := ψ.toRingHom.toAlgebra;
    <diamond 等式>` の形で直接書こうとしたところ、`differentIdeal V1
    Wn1` が要求する `Module.IsTorsionFree V1 Wn1`(`ψ` の単射性からしか
    出ない、`ψ` 選択**後**にしか証明できない事実)の instance 検索が、
    `∃`/`letI` を含む**シグネチャの型検査自体**(証明本体に入る前)で
    走ってしまい、`maxHeartbeats 2000000`(168秒)でも終わらないことを
    実測した。`hsep` で使った「`Σ'` で instance ごと束ねる」策も、
    今回は `Module.IsTorsionFree V1 Wn1` が `ψ`(存在文の中の束縛変数)
    に依存するため使えない(`ψ` を `Prop` の `∃` から取り出す
    `obtain` は目標が `Type` だと `Exists.casesOn` の
    `Prop`-限定に阻まれる、というさらに別の障害も確認した)。
    **結論: 戻り値は `True` のまま(内部で `hdiamond` という名の
    genuine な項を構築し、その型を docstring に逐語で記録する)方式
    に倒した**——数学的内容は完全に確立しているが、「他の証明から
    直接呼び出せる補題」としての再パッケージングは未解決のまま次回に
    持ち越す(方針候補: `ψ` の単射性を含む**全ての ψ 依存 instance を
    シグネチャの外側の明示的仮定として渡してもらう**——
    `differentIdeal_tower_diamond` 自身がそうしているのと同じ形——か、
    非 instance の明示引数として `@` で全て手渡しする、のいずれか)。

    ★これで Theorem 1.2 の証明骨格(3a: monogenicity・3b: 長さの評価・
    3c: `V_n` 塔の構成)のうち、**「different の塔の公式」という中核
    道具は完全に具体化できた**。残るのは (i) `cancel_conductor_delta`
    の `hspan_eq` をこの具体的 diamond 等式から導く接続、(ii) 再帰的な
    `V_n` 族(`V_1, V_2, ...`)を実際に構成すること——このセッションで
    確立した `V0→V1` の1ステップを繰り返し適用する形になるはずだが、
    「πの選び方(次の一様化元)」の再帰的な選定は未着手。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**(i) `hspan_eq`
    接続に向けて3つの補題を確立した**(`falt1_adjoin_top_of_finrank_eq`・
    `falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly`・
    `falt1BaseChangeGeneratorFull`、`KaehlerAux.lean`、commit
    `53152c9d`、`lake build`・`#print axioms` 確認済み・sorry 無し)。
    道筋の確認(数学的には閉じている):
    - `hspan_eq : spanDeriv = Idiff` は、両辺とも
      `Ideal.span {n · x^(n-1)}`(`x` は `Wₙ₊₁` の生成元、`ψ w = x`)
      に簡約されることを確認した——`Idiff = Ideal.map ψ (differentIdeal
      V0 V1) = Ideal.map ψ (span{n·w^(n-1)})`(`differentIdeal_eq_
      span_derivative` を `w` に適用)`= span{n·ψ(w)^(n-1)} =
      span{n·x^(n-1)}`(環準同型なので `n·(-)^(n-1)` を素通しする)、
      `spanDeriv = span{aeval x (deriv(minpoly Wₙ x))} =
      span{n·x^(n-1)}`(`conductor_mul_differentIdeal` を `x` に適用、
      `minpoly Wₙ x = X^n-Cπ'` なので微分は `n·X^{n-1}`)——**両辺が
      文字通り同じ式に簡約される**ので `rfl` に近いはずだと判明した。
    - `differentIdeal_eq_span_derivative`/`conductor_mul_differentIdeal`
      を実際に**呼び出す**には、`w`(`falt1BaseChangeGeneratorFull` の
      生成元)が field-level(`FractionRing V0` 上)でも全体を生成する
      ことを要求する——ring-level の `Algebra.adjoin V0 {w} = ⊤`
      だけでは足りない(これが `hspan_eq` 接続の技術的な核心だった)。
      これを `falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly` で解決:
      `minpoly.isIntegrallyClosed_eq_field_fractions'` で ring-level
      `minpoly V0 w = X^n-Cπ` を field-level `minpoly (FractionRing V0)
      (w:L) = fK` に持ち上げ、`falt1_adjoin_top_of_finrank_eq`(次元
      カウント: `Algebra.adjoin.powerBasis'` の次元 `=
      (minpoly K x).natDegree` と `Submodule.eq_top_of_finrank_eq`)
      で `Algebra.adjoin (FractionRing V0) {(w:L)} = ⊤` を得た。
    - `w`・`x` が「同じ `e1`・`e2` 由来」であることを保証するため
      (別々に `falt1GeneratorPackage` を呼ぶと `obtain` した witness
      の一致が保証されない)、`falt1BaseChangeGeneratorFull` で `ψ`・
      `w`・`x`・両方の生成元性を1つの証明にまとめた。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**上記の
    instance mismatch を解決し、`hspan_eq` を完全に確立した**
    (`falt1_hspan_eq`、commit `ea63551e`、`lake build`・
    `#print axioms` 確認済み・sorry 無し)。**真因は診断が違っていた**
    ——`@`明示引数や`show`強制ではなく、`differentIdeal_eq_span_
    derivative`/`conductor_mul_differentIdeal`を`AdjoinRoot fK`
    (条件付きで体になる型)経由の`w`・`x`に呼び出す**呼び出し側の
    文脈**に、`Fact (Irreducible fK)`・`FiniteDimensional
    (FractionRing V0)(AdjoinRoot fK)`・`Algebra.IsSeparable
    (FractionRing V0)(AdjoinRoot fK)`が**存在しなかった**ことが原因
    だった——これらは`falt1BaseChangeGeneratorFull`等のヘルパー内部で
    `haveI`されるだけで外に漏れないため、呼び出し側で**同じ3つを
    再構築**するだけで(`isDefEq`timeoutという紛らわしい症状のわりに)
    一瞬(3秒)で解決した。
    `conductor_mul_differentIdeal`の右辺`span{aeval x (deriv(minpoly
    Wₙ x))}`と`differentIdeal_tower_diamond`の`Idiff`
    (`Ideal.map ψ (differentIdeal V0 V1)`)が、両方とも
    `span{n·x^(n-1)}`に簡約されることを示した——
    `differentIdeal_eq_span_derivative`を`w`に適用して`differentIdeal
    V0 V1 = span{n·w^(n-1)}`を得て、環準同型`ψ`でmapし`ψ(w)=x`を使う
    だけ。**これでTheorem 1.2 item 3cの中核技術的障害は解決した。**
    残る接続: (a) `Idiff≠0`(`differentIdeal_ne_bot`——
    `[Module.Finite A B][Algebra.IsSeparable (FractionRing A)
    (FractionRing B)]`から`differentIdeal A B≠⊥`、`Ideal.map`の単射性
    保存と合わせれば出せる見込み、`Algebra.IsSeparable (FractionRing
    V0)(FractionRing V1)`という新しい instance が必要で未確立)を足して
    `cancel_conductor_delta`を実際に呼び出すこと、(b) 再帰的`V_n`族の
    構成。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**(a)を
    実行し、`cancel_conductor_delta`を実際に呼び出した**
    (`falt1_cancelConductorDelta_assembled`・`falt1_hsepV0V1_bundled`、
    commit `50b8564d`、`lake build`・`#print axioms`確認済み・sorry
    無し)。`Algebra.IsSeparable (FractionRing V0)(FractionRing V1)`
    は`falt1_hsep_bundled`の`V0→V1`単独版(`Wn`を経由しないぶん
    遥かに単純)として構築した。`differentIdeal_ne_bot`の呼び出しで
    もう1つ diamond の再来に遭遇した——`differentIdeal_tower_diamond`
    の`hsep`同様、`Algebra.IsSeparable`の下敷きとなる`Algebra`instance
    が`FractionRing.liftAlgebra`限定で、ambientな bracket instance
    では**instance検索自体が失敗する**(型mismatchですらなく
    「failed to synthesize」という紛らわしい形)——`@`明示引数で
    正しいinstanceを直接渡せば一瞬で解決した。
    **これで得られた結論**:
    ```
    conductor Wₙ x * differentIdeal V1 Wₙ₊₁
      = Ideal.map (algebraMap Wₙ Wₙ₊₁) (differentIdeal V0 Wₙ)
    ```
    (`Jₙ`が消去された、`δₙ`・`δₙ₊₁`を結ぶ関係式)——**これが
    Theorem 1.2・3b/3cの中核の代数的関係式であり、item 3cの技術的な
    障害はこれですべて解決した**。残るのは (i) 上式を`delta_tendsto_
    zero`が要求する`hrec`(長さの不等式)へ変換する
    ——`Module.length`の計算(3b後半、conductorの長さ評価)と
    `conductor(Wₙ,x)`が単項イデアルであることの利用が必要——、
    (ii) 再帰的`V_n`族(`V_1, V_2, ...`)を実際に構成すること
    (「πの選び方(次の一様化元)」の再帰的な選定は未着手)。

    ★(i)の下調べ(2026-09-04、未着手のまま次回に持ち越し):
    `length_quotient_span_singleton_mul`(既存、単項イデアルの場合の
    長さの加法性)を`conductor Wₙ x * differentIdeal V1 Wₙ₊₁`に適用
    するには`conductor Wₙ x`(`Ideal Wₙ₊₁`)が**単項**である必要がある
    ——これは`Wₙ₊₁`が**PID(局所化すれば DVR)**であれば任意のイデアルが
    単項なので自動的に満たされる。`Wₙ₊₁ := integralClosure Wₙ
    (AdjoinRoot gK)`は`gK`が Eisenstein なので古典的には**全分岐
    (totally ramified)**——`Wₙ`上唯一の素イデアルしか持たず、
    Dedekind整域+局所 ⟹ DVR、のはず——だが、これを mathlib で示す
    経路は未調査。★このファイル既存の教訓(「3c: 戦略転換」の節)が
    まさにこの種の「Eisenstein ⟹ 整閉・全分岐」を直接示そうとして
    「局所化・付値延長という mathlib にまだ薄い領域」に何度も
    突き当たったと記録している——**同じ壁が(i)でも予想される**。
    次回はまず`IsDiscreteValuationRing (integralClosure Wₙ
    (AdjoinRoot gK))`(または同値な「唯一の極大イデアル」)が
    mathlibの既存資産だけで示せるか、`node tools/decl-index.mjs
    --mathlib`で`ramificationIdx`・`IsDedekindDomain.isPrincipalIdealRing_
    of_...`・`IsLocalRing`周りを先に調査することから始めるのが
    効率的(証明を試す前に在庫を確認する、CLAUDE.mdの「在庫」原則)。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**予想は
    外れた——「同じ壁」には当たらなかった**(`falt1_
    isPrincipalIdealRing_of_finite_ext_of_DVR`、commit `a110a169`、
    `lake build`・`#print axioms`確認済み・sorry無し)。「全分岐
    (唯一の極大イデアル)」を**直接**示す必要は無く、**もっと弱い
    「極大イデアルが有限個」で十分**(`IsPrincipalIdealRing.of_
    finite_maximals`)と気づいたのが鍵だった。`Wₙ`が DVR(唯一の
    極大イデアル)であることと、有限拡大では上にある素イデアルが
    有限個であること(`IsDedekindDomain.primesOver_finite`)・整拡大
    では極大イデアルの引き戻しが極大であること(`Ideal.isMaximal_
    comap_of_isIntegral_of_isMaximal`)を組み合わせるだけで、`Wₙ₊₁`
    の極大イデアル全体が「有限集合の部分集合」だと分かり、部分集合の
    有限性だけで`IsPrincipalIdealRing Wₙ₊₁`が閉じた——「全分岐」という
    (証明がずっと重い)強い主張を経由しない近道が見つかったことになる。
    `Mathlib.RingTheory.DedekindDomain.PID`の import を追加した。
    次回はこの`IsPrincipalIdealRing`から`conductor(Wₙ,x)`の生成元を
    実際に取り出し(`IsPrincipalIdealRing`は存在の主張なので`obtain`で
    生成元を得る)、それが非零因子であること(`Wₙ₊₁`は整域なので
    非零⟹非零因子)を確認し、`length_quotient_span_singleton_mul`を
    実際に適用して長さの等式を得るところから続ける。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**実行し、
    長さの等式(`hlen_eq`)まで到達した**(`falt1_cancelConductorDelta_
    assembled`の拡張、commit `0bd39e3b`、`lake build`・`#print axioms`
    確認済み・sorry無し)。`conductor(Wₙ,x)≠0`は`length_quotient_
    span_singleton_mul`が要求する形とは別の道筋で出た——単項生成元を
    直接扱う代わりに、既に確立済みの`hcond`・`hspan_eq`・`hIdiff_ne`
    から`conductor(Wₙ,x)*differentIdeal Wₙ Wₙ₊₁ = Ideal.map ψ
    (differentIdeal V0 V1) ≠ 0`を作り、「積が非零なら整域では両因子
    とも非零」で`conductor(Wₙ,x)≠0`を導いた(`falt1_length_quotient_
    mul_of_ne_zero`が`I≠0`だけを仮定する形にしておいたので、単項
    生成元を明示的に扱う必要すら無かった)。得られた式:
    ```
    Module.length Wₙ₊₁ (Wₙ₊₁ ⧸ Ideal.map (algebraMap Wₙ Wₙ₊₁) (differentIdeal V0 Wₙ))
      = Module.length Wₙ₊₁ (Wₙ₊₁ ⧸ conductor Wₙ x) + Module.length Wₙ₊₁ (Wₙ₊₁ ⧸ differentIdeal V1 Wₙ₊₁)
    ```
    ★残る接続(次回持ち越し、これが真の核心の困難と評価する):
    この式を`hrec`(実数列の不等式)へ変換するには、
    (a) 左辺(`length(Wₙ₊₁/Ideal.map(differentIdeal V0 Wₙ))`)を
    `δₙ`(`Wₙ`側の量)と結ぶこと——`length_map_pow_of_ramificationIdx`
    (既存、素イデアルの冪の像の長さ=分岐指数倍)が使えるはずだが、
    `differentIdeal V0 Wₙ`を「`Wₙ`の極大イデアルの冪」の形に書き直す
    追加の手順が要る、
    (b) 右辺第2項(`length(Wₙ₊₁/differentIdeal V1 Wₙ₊₁)`)が
    `δₙ₊₁`そのもの、
    (c) 右辺第1項(`length(Wₙ₊₁/conductor(Wₙ,x))`)の**下界評価**——
    これが以前から「独立の古典的整数論(discriminant の塔・
    conductor-discriminant の関係)を要する」と評価してきた核心の
    困難に対応し、原文が「劣加水性・完全列」評価と呼ぶ箇所と一致する
    と見られる。まだ着手していない。

    ★2026-09-04、下調べで(a)の評価を訂正: `length_map_pow_of_
    ramificationIdx`の前提`Ideal.map(algebraMap Wₙ S)(maximalIdeal
    Wₙ) = maximalIdeal(S)^e`は`S`(ここでは`Wₙ₊₁`)が**単一の極大
    イデアルしか持たない(局所)**ことを暗黙に要求する——これは
    (i)の下調べで確立した`Wₙ₊₁`が**PID**であること(有限個の極大
    イデアル)より強い主張で、**「全分岐」の問題そのものに実質的に
    帰着する**と判明した。★★つまり(a)は(c)とは独立の「もう1つの
    核心の困難」ではなく、**同じ「全分岐 / discriminant の塔」の
    問題の別の顔**だった可能性が高い——`Wₙ₊₁`が実際に全分岐なら
    (a)(c)双方が(ラマヌジャン係数`e=n`と conductor-discriminant の
    関係から)同時に解決するはずで、逆に全分岐を回避する一般の
    (複数素点の)議論を(a)単独で組む必要は無いと見てよい。
    **結論**: 次回はまず「`Wₙ₊₁`が全分岐である」こと自体
    (`(maximalIdeal Wₙ).primesOver Wₙ₊₁`が単集合であること)を
    正面から示すことに一本化するのが筋が良い——`Ideal.sum_
    ramification_inertia_eq_finrank`(既存、`Σ e_i·f_i = [Wₙ₊₁:Wₙ]`)
    と、Eisenstein多項式の次数`n`=`[Wₙ₊₁:Wₙ]`を組み合わせれば、
    「剰余体拡大次数`f=1`かつ素点が1つ」を`n = e·f`の等式の一意分解
    (Wₙの剰余体が有限とは限らないので単純な数え上げでなく、
    Eisenstein性から直接`e≥n`を出す議論が必要になる見込み)から
    追い込める可能性がある——これが3c本来の「全分岐」証明そのもの
    であり、当初「3c: 戦略転換」で回避した論点に結局戻ることになる。

    ★2026-09-04、近道の探索を完了(見つからず、評価を確定):
    `Module.length`の base change 公式(`IsLocalRing.length_
    baseChange`・`CovBy.length_baseChange`、`RingTheory/LocalRing/
    Length.lean`)は一見「全分岐を回避できるかもしれない」候補
    だったが、**どちらも `[IsLocalRing B]`(拡大先の局所性)を要求
    する**と`#check`で確認した(mathlib-index.txt の抜粋だけでは
    この仮定が見えず、一時的に「回避できるかも」と誤認した——
    `#check`で完全な型を見ることの重要性の再確認)。CRT で複数の
    局所成分に分解して長さを足し合わせる一般補題も mathlib に
    見当たらなかった。**結論(確定)**: (a)(c)の接続は、`Wₙ₊₁`が
    全分岐であることの直接証明を回避できず、これが Theorem 1.2の
    残り(3b後半・3c)における唯一かつ真の核心の困難である
    ——これ以上の近道の探索は費用対効果が低いと判断し、次回は
    (探索ではなく)全分岐の証明そのものに着手するか、この論点を
    独立の未解決事項として記録した上で3cの別の側面(再帰的`V_n`族の
    構成)に進むか、を判断してから着手する。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**(c)は
    実は解消していた——重要な発見(要検討の「逸脱」)**(commit
    `4de9d57c`)。`hxadjoin`(`x` が `Wₙ₊₁` を `Wₙ`-**代数として**
    生成する、`Wₙ[x]=Wₙ₊₁`、既に確立済み)から`conductor_eq_top_
    iff_adjoin_eq_top`で`conductor(Wₙ,x)=⊤`が直ちに従い、
    `length(Wₙ₊₁/conductor(Wₙ,x))=0`——(c)で懸念していた「下界評価」
    は不要で、`hlen_eq`は**補正項の無い完全な等式**になる:
    ```
    length(Wₙ₊₁/Ideal.map(algebraMap Wₙ Wₙ₊₁)(differentIdeal V0 Wₙ))
      = length(Wₙ₊₁/differentIdeal V1 Wₙ₊₁)
    ```
    ★★★ただし**これは「このセッションの構成固有の単純化」である
    可能性が高いと評価する**——`Wₙ₊₁` を `V1` と**全く同じ**
    Eisenstein多項式`gK`の base change から作る、という選択そのもの
    が、`x` の ring-level 生成性(Eisenstein性から monogenicity が
    「タダで手に入る」、`adjoin_eq_top_in_integralClosure_of_
    isEisensteinAt`)を保証してしまい、conductor が自明になった
    ——Faltings原論文の Theorem 1.2 の一般論では `Wₙ` は(その前段の
    再帰で構成された)**任意の**代数拡大でよく、`x` が必ずしも
    `Wₙ₊₁` をring-levelで生成する保証は無い可能性がある(そここそが
    原文の「非正規性の評価(劣加水性・完全列)」が本来扱う核心の
    困難だったはず)。**次回の要検討事項(「逸脱」の候補、確定次第
    `ResearchPaper/`か本ファイルに正式記録すること)**: 再帰的`V_n`
    塔を実際に構成する際、`Wₙ`(=前段で構成された`V_{n}`の役割)が
    本当に「`V_{n+1}`と同じEisenstein多項式のbase change」という
    形で現れるのか、それとも一般にはそうでないのかを、原文
    (Theorem 1.2の証明、`falt1-goal.md`冒頭の物理ページ情報参照)
    に立ち戻って確認すること——もし原文の設定で本当に`Wₙ`がこの形
    (V_{n+1}の base change)なら、この単純化は「逸脱」ではなく
    **原文に忠実な特殊化**であり、conductor が消えるのは原文の
    数学的事実そのものということになる(要検証)。

    ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、**原文に
    立ち戻って確認した——重要な訂正**。`Faltings - p-adic Hodge
    theory.txt`(物理p.257=印字p.223付近、Chapter I Theorem 1.2の
    原文)と`Brinon Conrad - CMI Summer School Notes...txt`
    Exercise 13.7.4(物理p.230付近)を読み直した結果:
    - `Wₙ`(原文の`B_n`)は**固定された1つの拡大`W`(原文の`B`)**を
      `V_n`(原文の`A_n`)に沿って base change して正規化したもの
      (`W_n := normalization(V_n ⊗_V W)`)——**`V_{n+1}`と同じ
      Eisenstein多項式の base change ではない**。私が「`Wₙ`」と
      呼んできたものは、この`W_n`の(半局所環の)**1つの成分**
      (原文が明記する「or more precisely of each factor of the
      semilocal ring」)に相当する。
    - ただし`V_{n+1}=V_n[x]/(f)`(`f`がEisenstein)なら、
      `W_n ⊗_{V_n} V_{n+1} = W_n[X]/(f)`(**同じ`f`**をbase change
      したもの)となり、正規化前の構成は**私の`Wₙ₊₁ := integralClosure
      Wₙ(AdjoinRoot gK)`と完全に一致する**——これは「逸脱」ではなく
      **原文に忠実**だったと確認できた。
    - **重要な訂正点**: 原文自身が「`W_n`は`A_n`上有限な離散付値環
      **有限個の積**(半局所環)」と明記しており、`Wₙ₊₁`が全分岐
      (単一素点)であることは**原文でも一般には主張されていない**
      ——むしろ**複数の素点に分裂しうることが前提**になっている。
      これは、私が(a)の接続のために「全分岐」を証明しようとしていた
      方向性が**そもそも誤りだった可能性が高い**ことを意味する
      ——mathlibで見つからなかったのは「まだ薄い」からではなく、
      **一般には成り立たない主張だから探しても見つからない**
      可能性がある。PID性(有限個の極大イデアル)の証明で「全分岐を
      回避できた」ことは、実は**原文の設定(半局所環)に合致した
      正しい一般化のレベル**だったと分かる。
    - 原文自身の議論(Exercise 13.7.4 (4)-(6))は「conductor」という
      語を使わず、**Kähler微分`Ω¹`の完全列**(`Ω¹_{Bn/An}⊗Bn+1 →
      Ω¹_{Bn+1/An} → Ω¹_{Bn+1/An+1}`)の kernel/cokernel の長さを
      直接評価する——このセッションで組み立てた`cancel_conductor_
      delta`ベースの経路は、原文と**別の(ただし関連した)証明経路**
      である可能性が高い。両者が本当に同値な結論を導くか、この
      `conductor=⊤`という発見が原文の議論とどう対応するかは、
      **次回、Exercise 13.7.4 (4)-(6)を精読してから再検討する**
      必要がある——特に、私の`Wn1`が半局所(複数成分)の場合に
      `conductor(Wn,x)=⊤`が各成分ごとに成り立つか(`x`が各成分を
      個別に生成するか)は未確認。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **Exercise 13.7.4 (4)の第一歩(完全列の左側の単射性)を確立
      した**(`falt1_mapBaseChange_injective_adjoinRoot`、commit
      `6c78b84c`、`lake build`・`#print axioms`確認済み・sorry無し)。
      Lemma 1.1 の土台で既に完成していた`mapBaseChange_injective_
      of_nzd`(`B=Polynomial V/(f)`かつ`f'`が非零因子の場合の
      `mapBaseChange`単射性)を、`V0→Wn→AdjoinRoot g`
      (`g:=X^n-Cπ'`、base changeする**前**のEisenstein多項式)に
      直接適用するだけで、
      ```
      KaehlerDifferential.mapBaseChange V0 Wn (AdjoinRoot g) :
        AdjoinRoot g ⊗_Wn Ω¹_{Wn/V0} → Ω¹_{AdjoinRoot g/V0}
      ```
      の単射性が得られた——これが Exercise 13.7.4 の完全列
      `0 → Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0} → Ω¹_{Wₙ₊₁/V0} → Ω¹_{Wₙ₊₁/Wₙ} → 0`
      の**左側**に相当する(`AdjoinRoot g ≃ₐ[Wₙ] Wₙ₊₁`という既存の
      同型`falt1AdjoinRootEquivIntegralClosure`を経由すれば`Wₙ₊₁`
      自身の単射性に持ち上がる見込み)。次回はこの transport
      (`omegaCongr`、既存)から続け、右側の全射性(標準的、一般の
      塔で成り立つはず)と合わせて完全列そのものを組み立て、
      Nakayama+elementary divisors(step 4後半)・discriminantの塔
      (step 5)へ進む。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **上記の transport を実際に完了し、`Wₙ₊₁`自身の単射性を確立
      した**(`falt1_mapBaseChange_injective_wn1`・
      `falt1_mapBaseChange_injective_wn1_concrete`、
      `Found/Falt1/KaehlerAux.lean`、commit `fe7520bd`、`lake build`
      完走・sorry無し)。`KaehlerDifferential.mapBaseChange_comp`・
      `map_comp_der`のような直接の自然性補題はmathlibに存在しない
      と確認済みだったが、自前でその自然性を証明し直す必要は無かった
      ——Lemma 1.1の土台で既に`falt1MapBaseChangeInjective`のために
      組み立てていた`mapBaseChange_injective_transport`(代数同型
      `e:A≃ₐ[V]B`に沿って`mapBaseChange`の単射性を輸送する、
      `omegaCongr`を内部で使う既存の補題)を、`A:=AdjoinRoot g`・
      `B:=Wₙ₊₁`・`e:=falt1AdjoinRootEquivIntegralClosure`として
      そのまま適用するだけで閉じた。「既存の資産を掘り直す前に
      在庫を確認すべきだった」という教訓——`tools/decl-index.mjs`
      でのgrepを`mapBaseChange_injective`という結論の語で先に
      引いていれば1手で見つかっていた。
      これでExercise 13.7.4 step (4)の完全列の**左側**(単射性)が
      Falt1の実際の`Wₙ₊₁`に対して完成した。次回: 右側の全射性
      (`KaehlerDifferential.mapBaseChange_surjective`等、標準的な
      事実のはず・未着手)を組み立てて完全列そのものを確立し、
      Nakayama+elementary divisors(step 4後半のkernel評価)・
      discriminantの塔(step 5のcokernel評価)へ進む。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **完全列の残り2成分が実は無条件でmathlibに既にあると判明し、
      長さの加法公式を確立した**(`falt1_kaehler_length_exact_wn1`、
      commit `388c961e`、`lake build`完走・sorry無し)。第二基本完全列
      `0 → Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0} → Ω¹_{Wₙ₊₁/V0} → Ω¹_{Wₙ₊₁/Wₙ} → 0`
      の中央の完全性(`KaehlerDifferential.exact_mapBaseChange_map`)・
      右側の全射性(`KaehlerDifferential.map_surjective`、`map R S B B`
      の形——`Ω[B/R]→Ω[B/S]`が**任意の`R→S→B`の塔で無条件に全射**)は
      分岐や分離性の仮定なしに成り立つ一般論としてmathlibに既にあった。
      左側の単射性(直前のマイルストーン)と`Module.length_eq_add_of_exact`
      を合わせるだけで
      ```
      length_{Wₙ₊₁}(Ω¹_{Wₙ₊₁/V0}) =
        length_{Wₙ₊₁}(Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0}) + length_{Wₙ₊₁}(Ω¹_{Wₙ₊₁/Wₙ})
      ```
      が閉じた。これがExercise 13.7.4 step (4)-(6)の議論の土台
      (kernel・cokernelの長さ評価に相当する右辺2項)——次回はこの
      右辺2項をそれぞれ評価する: 第1項(kernel)はstep (4)後半の
      Nakayama+elementary divisors(`Ω¹_{Wₙ/V0}`が`d+1`個の生成元を
      持つことと組み合わせる)、第2項(`Ω¹_{Wₙ₊₁/Wₙ}`、cokernelに
      相当)はstep (5)のdiscriminantの塔。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **右辺第2項を実際に`differentIdeal`の言葉に変換した**
      (`falt1_kaehler_length_exact_wn1_cokernel`、commit `84793b74`、
      `lake build`完走・sorry無し)。Lemma 1.1(`falt1CokernelLengthEq`、
      `Ω_{W/V}`の長さ`=W⧸differentIdeal V W`の長さ)を`V:=Wₙ`・
      `W:=Wₙ₊₁`に適用するだけで
      ```
      length_{Wₙ₊₁}(Ω¹_{Wₙ₊₁/V0}) =
        length_{Wₙ₊₁}(Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0}) +
        length_{Wₙ₊₁}(Wₙ₊₁ ⧸ differentIdeal Wₙ Wₙ₊₁)
      ```
      **重要な発見**: この右辺第2項`length(Wₙ₊₁⧸differentIdeal Wₙ Wₙ₊₁)`
      は、`cancel_conductor_delta`経由で既に確立していた`hlen_eq`
      (`falt1_cancelConductorDelta_assembled`)にも**全く同じ量**
      `differentIdeal V1 Wn1`(`≠ Wₙ Wₙ₊₁`だが同一のもの、`V1 ≃ₐ[V0] Wₙ`
      の下で表現が違うだけ——要再確認)として現れる。つまり
      **2つの証明経路(conductor経由・Kähler微分の完全列経由)は
      独立ではなく、同じ`differentIdeal`の長さを介して繋がっている**
      ことが確認できた——両者が同値かは未確定だが無関係ではない。
      残るは右辺第1項`length(Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0})`(kernel側)の
      評価のみ——次回はこれに着手する。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **右辺第1項(kernel側)も評価し終えた**
      (`falt1_kaehler_length_exact_wn1_kernel`、commit `0292d557`、
      `lake build`完走・sorry無し)。Lemma 1.1(`falt1CokernelIsoLinear`、
      `Ω¹_{Wₙ/V0} ≃ₗ[Wₙ] Wₙ⧸differentIdeal V0 Wₙ`)を`Wₙ→Wₙ₊₁`へ
      `LinearEquiv.baseChange`するだけで
      `Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0} ≃ₗ[Wₙ₊₁] Wₙ₊₁⊗_{Wₙ}(Wₙ⧸differentIdeal V0 Wₙ)`
      が出るが、これをさらに`Wₙ₊₁⧸(differentIdeal V0 Wₙ).map(...)`に
      単純化する一般補題`A⊗[R](R⧸I)≃ₐ[A]A⧸I.map(algebraMap)`が
      mathlibに直接無く(`Algebra.TensorProduct.tensorQuotientEquiv`
      と`rid`を自分で合成する必要があった)、新規に
      `falt1_tensorQuotientEquiv_algebraMap`として確立した。結果:
      ```
      length_{Wₙ₊₁}(Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0}) =
        length_{Wₙ₊₁}(Wₙ₊₁ ⧸ (differentIdeal V0 Wₙ).map(algebraMap Wₙ Wₙ₊₁))
      ```
      **これで`falt1_kaehler_length_exact_wn1`の右辺2項がどちらも
      `differentIdeal`の言葉に変換できた**:
      ```
      length_{Wₙ₊₁}(Ω¹_{Wₙ₊₁/V0}) =
        length_{Wₙ₊₁}(Wₙ₊₁⧸(differentIdeal V0 Wₙ).map(algebraMap Wₙ Wₙ₊₁)) +
        length_{Wₙ₊₁}(Wₙ₊₁⧸differentIdeal Wₙ Wₙ₊₁)
      ```
      **重要**: 右辺第1項の`Wₙ₊₁⧸(differentIdeal V0 Wₙ).map(algebraMap
      Wₙ Wₙ₊₁)`は`cancel_conductor_delta`経由で既に確立していた
      `hlen_eq`(`falt1_cancelConductorDelta_assembled`、簡約後:
      `length(Wₙ₊₁⧸Ideal.map(algebraMap Wₙ Wₙ₊₁)(differentIdeal V0 Wₙ))
      = length(Wₙ₊₁⧸differentIdeal V1 Wₙ₊₁)`)の**左辺そのもの**——
      2つの証明経路が完全に噛み合う点まで来た。次回はこれを実際に
      代入し、
      ```
      length_{Wₙ₊₁}(Ω¹_{Wₙ₊₁/V0}) =
        length_{Wₙ₊₁}(Wₙ₊₁⧸differentIdeal V1 Wₙ₊₁) +
        length_{Wₙ₊₁}(Wₙ₊₁⧸differentIdeal Wₙ Wₙ₊₁)
      ```
      という**完全に`differentIdeal`だけで書かれた閉じた式**を導出する
      ——これがTheorem 1.2の再帰(`δ_n→0`)の核心になりうる。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      実際に`falt1_cancelConductorDelta_assembled`の結論を`True`から
      この最終式に置き換えようとしたが、**着手して初めて分かった
      追加の障害**を記録する: `falt1_kaehler_length_exact_wn1_kernel`/
      `_cokernel`(Lemma 1.1 の`falt1CokernelIsoLinear`/`LengthEq`を
      経由)を呼ぶには、`[IsIntegralClosure Wn V0 (FractionRing Wn)]`
      (`Wn`が`V0`の`FractionRing Wn`内での整閉包であること)・
      `[IsIntegralClosure Wn1 Wn (AdjoinRoot gK)]`という**この
      assembled定理では まだ確立していないインスタンス**が要る
      ——`Wn`は一般に与えられた`IsDiscreteValuationRing`(それ自体が
      積分閉整域)で`Module.Finite V0 Wn`も仮定にあるので、一般論
      (整閉整域+有限次数拡大⟹整閉包)から出るはずだが、mathlibでの
      具体的な組み立て方(`IsIntegralClosure.of_isIntegrallyClosed`
      のような名前の補題を探すところから)は未確認。**「軽い代入で
      閉じる」という見立ては誤りで、独立の一歩(この整閉包インスタンス
      の確立)が要る**——このセッションで繰り返し見たパターン
      (`length(R/(I·J))=...`の加法性の時と同型)がここでも再現した。
      安全のため`falt1_cancelConductorDelta_assembled`の型は`True`の
      まま(変更を revert 済み、`git diff`で無変更を確認)にして
      次回へ引き継ぐ——次回はまず
      `IsIntegralClosure Wn V0 (FractionRing Wn)`のインスタンスを
      独立に確立することから始めること。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      上記のインスタンス自体は**独立に解決できた**(REPL で確認済み、
      未コミット・調査のみ): `IsIntegralClosure.of_isIntegrallyClosed`
      (mathlib、`[IsIntegrallyClosed R][Algebra S R][Algebra S K]
      [IsScalarTower S R K][Algebra.IsIntegral S R]⟹IsIntegralClosure R S K`)
      を`R:=Wn, S:=V0, K:=FractionRing Wn`に適用するだけで
      `IsIntegralClosure Wn V0 (FractionRing Wn)`が(`hinjV0Wn`すら
      使わず)閉じた。`Wn1`側はさらに簡単で、`Wn1:=integralClosure Wn
      (AdjoinRoot gK)`という**構成そのもの**から`integralClosure.
      isIntegralClosure`(mathlib、常に成立)で無条件に`IsIntegralClosure
      Wn1 Wn (AdjoinRoot gK)`が出る。

      **しかし、着手して`falt1_kaehler_length_exact_wn1_kernel`を実際に
      `Wn`に適用しようとして、より深い構造的な障害を発見した**:
      Lemma 1.1(`falt1CokernelIsoLinear`/`LengthEq`)は`W`が`V`上
      **単項生成**(`w:W`で`hadjoin:V[w]=⊤`かつ`hw:K(w)=L`)であることを
      本質的に要求する。`Wn1`(`V1`上・`Wn`上どちらも)は`falt1BaseChange
      GeneratorFull`で生成元`x`が明示的に手に入るが、**`falt1_cancelConductorDelta_
      assembled`の中で`Wn`自体は`V0`上の抽象的な有限拡大として与えられて
      おり、単項生成元は用意されていない**(`Wn`は構成的に`Vₙ`として
      逐次構築されるはずだが、この定理はそこまで遡らず`Wn`を不透明な
      仮定として扱っている)。つまり kernel 項の評価には**`Wn/V0`の
      単項生成元(と、原始生成元の環レベル⟹体レベルの生成の遺伝、
      `falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly`のWn版に相当する
      未証明の一般事実)という、当初の想定より一段深い追加データ**が
      要ることが判明した——「軽い代入」という見立ては今回も誤りだった。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **上記の代替路を独立の一般補題として確立した**
      (`falt1_differentIdeal_tower_length`、commit `ce3743b0`、
      `lake build`完走・sorry無し)。`Wn`の単項生成元は一切不要——
      塔`A→B→C`に対して
      ```
      length_C(C/differentIdeal A C) = length_C(C/differentIdeal B C) +
        length_C(C/Ideal.map(algebraMap B C)(differentIdeal A B))
      ```
      が直接得られる、Kähler微分の完全列とは独立の経路。`A:=V0,B:=Wn,
      C:=Wn1`に適用し`falt1_cancelConductorDelta_assembled`の`hlen_eq`
      と組み合わせれば目的の式に到達できる見込みだが、**実際に
      `falt1_cancelConductorDelta_assembled`へ組み込もうとして別の
      障害に遭遇した**: この定理の結論(`differentIdeal A C`等)を
      `falt1_cancelConductorDelta_assembled`の**型シグネチャ**に書こう
      とすると、`IsDedekindDomain Wn1`・`Algebra V1 Wn1`等のインスタンス
      が必要だが、これらは(`Fact(Irreducible gK)`や`ψ`など)**証明本体
      の中でしか確立されない**ため、シグネチャの elaboration 時点では
      見つからない(`failed to synthesize instance`、実際に`lake build`
      で確認)。`Algebra V1 Wn1`は`ψ`(証明の中で`obtain`される)に依存
      するため、これをシグネチャで別途仮定すると`letI`で内部に再導入
      する`ψ`ベースの instance と**衝突しうる(diamond のリスク)**——
      このセッションで繰り返し見た「証明の中の局所知識をシグネチャに
      持ち上げようとして instance が壊れる」パターン。安全のため
      `falt1_cancelConductorDelta_assembled`自体の変更は revert し
      (`git diff`で無変更を確認)、`falt1_differentIdeal_tower_length`
      だけを独立の一般補題として commit した。次回: 最終結論は`V1`を
      経由しない形(`Ideal.map(algebraMap Wn Wn1)(differentIdeal V0 Wn)`
      のまま)で組み立て、`hlen_eq`で`V1`側に変換する一歩は
      `falt1_cancelConductorDelta_assembled`の**証明の中でだけ**行う
      よう設計し直すこと。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      上記の設計(`V1`を経由しない・`IsDedekindDomain`を追加の明示的
      仮定として渡す・`Ideal.top_mul`で`hIWnWn1_ne`を作る・
      `falt1_differentIdeal_tower_length hsep hIWnWn1_ne`で結ぶ)を
      実際に試したところ、**instance 探索ではなく elaboration の
      性能上の壁**に当たった: `set_option maxHeartbeats 3000000`
      (既にこの定理に設定済み、通常の150倍)まで上げても
      `(deterministic) timeout at isDefEq` で失敗する(`lake build`
      実測、479秒かけてタイムアウト)。原因と見立てられるのは、
      `Wn1`(`integralClosure Wn (AdjoinRoot(...))`という巨大な入れ子
      式)を**シグネチャに`set`なしでそのまま3回以上書いた**ため、
      `differentIdeal`・`Ideal.map`の型検査のたびに同じ巨大な式を
      unfold して instance 探索し直す必要が生じ、証明本体側の
      `set Wn1 := ...`(略記)と噛み合わないこと。安全のため
      `falt1_cancelConductorDelta_assembled`への変更は再度 revert し
      (`git diff`で無変更を確認)、`falt1_differentIdeal_tower_length`
      (独立の一般補題、既に commit 済み)はそのまま活かした。
      **次回の具体的な方針**: シグネチャレベルで`Wn1`を`let`束縛
      (`theorem foo ... : let Wn1 := ...; <本文>`という形)して1回だけ
      定義し、`intro`で導入してから使うか、あるいは`falt1_
      cancelConductorDelta_assembled`自体を分割し、`Wn1`・`V1`を
      **具体的な`integralClosure`式ではなく抽象的な型変数として渡す
      新しいラッパー定理**(`falt1_differentIdeal_tower_length`と
      同じ抽象化レベル)を用意して、そちらに`hlen_eq`等の**証明済み
      の事実だけ**を引数として渡す設計にすること——巨大な具象型を
      シグネチャに繰り返し書く設計そのものを避けるのが本筋。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      この方針を実行し、`falt1_differentIdeal_tower_length_via_hlen_eq`
      (commit `c4006562`、`lake build`完走・sorry無し)として確立した
      ——`falt1_differentIdeal_tower_length`を、「`Ideal.map(...)`の項が
      既に別のイデアル`E`の長さと一致することが分かっている」場合に
      特化した3行のコロラリー(`hlen_eq`を渡すだけ)。抽象的な型変数
      `A,B,C`のまま(具象の`integralClosure`式を経由しない)なので、
      前回のelaboration timeout(479秒)の原因(巨大な式をシグネチャに
      繰り返し書く)を回避——今回は数秒で`lake build`が通った。
      **残る作業**: これを実際に`falt1_cancelConductorDelta_assembled`
      の**証明の中で**(シグネチャではなく)`A:=V0,B:=Wn,C:=Wn1,
      E:=differentIdeal V1 Wn1`として適用するのは可能なはず(証明の
      中なら`Wn1`は`set`で略記済みで高速)。ただし`falt1_
      cancelConductorDelta_assembled`自体の**戻り値の型**(現在
      `True`)を意味のある式に変える際は、今回と同じelaboration
      timeoutの罠に注意すること——`V0,Wn`だけを引数に取り、戻り値の
      型で`Wn1`に触れる必要がある場合は`∃ Wn1, ...`のような存在型に
      するか、証明の最初に`obtain`で`Wn1`を取り出してから`show`で
      書き直す、といった設計が必要——次回はここから続ける。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **上記の設計を実際に成功させ、`differentIdeal`だけで書かれた
      閉じた式を完成させた**(`falt1_theorem12_differentIdeal_length`、
      commit `85a2e5fa`、`lake build`完走(224秒)・`#print axioms`確認・
      sorry無し)。timeoutの根本原因は「`Wn1`(巨大な入れ子式)を
      シグネチャに`set`なしで複数回ベタ書きすると、`differentIdeal`・
      `Ideal.map`の型検査のたびにinstance探索をやり直す」ことだった
      ——これを**名前付き`def`**(`Falt1Wn1 V0 Wn π n`)で解決した:
      1. `Falt1Wn1`を1回だけ`def`として定義。
      2. 必要な基本instance(`CommRing`・`Algebra Wn`・`Algebra V0`・
         `IsScalarTower V0 Wn`)を`inferInstanceAs`で1回だけ登録
         (`Falt1Wn1`は`noncomputable def`=半簡約なので自動では
         unfoldされない、これが鍵)。
      3. 定理のシグネチャで`Falt1Wn1 V0 Wn π n`を(何度書いても)
         短い適用として型検査——実測で数秒(以前は479秒timeout)。
      4. 証明の中では通常通り`set Wn1 := integralClosure Wn
         (AdjoinRoot gK)`で局所展開し、`Falt1Wn1`側のinstance仮定は
         `‹...›`(assumption)で`Wn1`側へ1回だけ変換。
      5. 最後は`have hresult : (Wn1で書いた具体的な型) := ...`と
         明示的に型注釈してから`exact hresult`——`exact`の型推論を
         `hsep`(Wn1型を確定させる)に**先に**決めさせることで、
         ゴール(`Falt1Wn1`型)との defeq 判定が(既に確定した局所
         instanceの再利用だけで済むので)高速になる。`show`(ゴールを
         生の式に書き換える)を使わない・型注釈付き`have`で先に確定
         させる、という2点が今回の鍵だった。

      得られた結果(`falt1_cancelConductorDelta_assembled`と同じ
      導出を`falt1_differentIdeal_tower_length`で結んだだけ):
      ```
      length_{Wₙ₊₁}(Wₙ₊₁⧸differentIdeal V0 Wₙ₊₁) =
        length_{Wₙ₊₁}(Wₙ₊₁⧸differentIdeal Wₙ Wₙ₊₁) +
        length_{Wₙ₊₁}(Wₙ₊₁⧸Ideal.map(algebraMap Wₙ Wₙ₊₁)(differentIdeal V0 Wₙ))
      ```
      **これがTheorem 1.2・3b/3cの中核の代数的関係式の、`differentIdeal`
      だけで書かれた完全な閉じた形**——`falt1_cancelConductorDelta_
      assembled`が`True`しか返せなかった制約をついに解消した。

      ★次回の残る作業: 左辺`differentIdeal V0 Wₙ₊₁`を`Ω¹_{Wₙ₊₁/V0}`の
      長さ(Lemma 1.1、`falt1CokernelLengthEq`経由)と結びつけるには
      `Wₙ₊₁`の`V0`上の単項生成元(体レベルでも生成する元)が要る——
      これは以前発見した障害(`falt1_kaehler_length_exact_wn1_kernel`
      をWnに適用しようとして遭遇したのと同じ種類の障害だが、今回は
      `Wₙ₊₁`自身についてなので、実は`falt1BaseChangeGeneratorFull`の
      `x`(`Wₙ`上の生成元)とは別に、`V0`上の生成元を別途探す必要が
      ある)。あるいは、`Ω¹_{Wₙ₊₁/V0}`を経由せず`differentIdeal V0
      Wₙ₊₁`という量そのものを直接評価する経路(discriminantの塔、
      未着手)に切り替える方が近道の可能性もある——次回はどちらが
      有望か検討することから始める。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **単項生成元を経由しない代替路を実際に見つけ、確立した**
      (`falt1_differentIdeal_diamond_length`、commit `408f3081`、
      `lake build`完走・sorry無し)。`differentIdeal V0 Wₙ₊₁`を`Ω¹`の
      長さに結びつけようとして単項生成元の壁にぶつかったが、そもそも
      **その必要が無い**ことに気づいた: `differentIdeal_tower_diamond`
      は`differentIdeal V0 Wₙ₊₁`を`V1`経由・`Wₙ`経由の**2通り**で
      計算した結果が一致するという主張(`Vₙ₁`・`Wₙ`側の単項生成元は
      どちらも既に不要)なので、その**両辺**に`falt1_length_quotient_
      mul_of_ne_zero`を適用するだけで、`Ω¹_{Wₙ₊₁/V0}`を一切経由せず
      直接
      ```
      length_{Wₙ₊₁}(Wₙ₊₁/differentIdeal V1 Wₙ₊₁) +
        length_{Wₙ₊₁}(Wₙ₊₁/Ideal.map(algebraMap V1 Wₙ₊₁)(differentIdeal V0 V1))
        =
      length_{Wₙ₊₁}(Wₙ₊₁/differentIdeal Wₙ Wₙ₊₁) +
        length_{Wₙ₊₁}(Wₙ₊₁/Ideal.map(algebraMap Wₙ Wₙ₊₁)(differentIdeal V0 Wₙ))
      ```
      という等式が出た。これは`Jₙ:=differentIdeal Wₙ Wₙ₊₁`(現在の段の
      非正規性)と`Jₙ₊₁:=differentIdeal V1 Wₙ₊₁`(次の段の非正規性、
      このセッションの単段indexingでは`V1`が`Vₙ₊₁`役)を結ぶ、
      **段をまたぐ長さの等式**——Theorem 1.2の再帰(`δₙ→0`)に必要と
      見込まれる関係式そのものである可能性が高い。

      **設計上の判断**: この定理は`Vₙ,Vₙ₁,Wₙ,Wₙ₁`すべてを抽象的な
      型変数のまま(`falt1_differentIdeal_tower_length`と同じ抽象化
      レベルで)確立した。`falt1_theorem12_differentIdeal_length`の
      ように具体的な`Falt1V1`/`Falt1Wn1`へ具体化した「閉じた」版を
      作ろうとすると、`Algebra V1 Wn1`(`ψ`という非自明な選択に依存)を
      シグネチャに固定する必要があり、これは`ψ`が証明本体でしか
      手に入らない(多くの補助仮定を要する)ため、`falt1_theorem12_
      differentIdeal_length`で最初に遭遇したのと同じ instance 衝突の
      罠に落ちる。**抽象的なままにしておくのが正しい設計**だと判断
      した——将来、再帰的な`V_n`塔を実際に構成するときには、各段が
      自前の自然な`Algebra Vₙ₁ Wₙ₁`インスタンスを持つはずなので、
      この一般補題をそのまま各段に適用できる。次回: この関係式を
      実際に「`δₙ`が減少する」という主張(3bの実数列の不等式、既に
      証明済みの`hrec`前提部分)へどう繋げるか、また`differentIdeal
      Wₙ Wₙ₊₁`(`Jₙ`)自体の評価(step 5、discriminantの塔)へ進む
      ことを検討する。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **`hrec`への接続を実際に調べ、残る障害を2つとも具体的に特定した**
      (`delta_tendsto_zero`(line 1518近辺、完成済み)の`hrec : ∀n, δ(n+1)
      ≤ δn - min(1,δn/(d+1))/(d+2)`へ接続する経路を検討):
      1. **`δn`を実数として定義する経路そのものが2通りあり、どちらも
         未完成の別の道具を要する**:
         - (a) `δn := Module.length(...)`を`Ω¹`経由で定義する場合、
           Lemma 1.1(`falt1CokernelLengthEq`)で`Wₙ₊₁`の`V0`上の
           単項生成元が要る(既出の障害)。
         - (b) `δn`を「`Wₙ₊₁`の極大イデアルの冪としての`differentIdeal`
           の**指数**」(Faltings原文の記法)として定義する場合、
           `length_map_pow_of_ramificationIdx`(完成済み、line 1775
           近辺)が使えるが、これは`Wₙ₊₁`が**単一のDVR**(全分岐)
           であることを要求する。mathlibで
           `Polynomial.IsEisensteinAt`(既存)と「局所性・全分岐」を
           結ぶ補題を探したが**見当たらなかった**(`grep`で確認、
           `IsEisensteinAt`関連の補題はirreducibility・adjoin関連の
           みで、Henselの補題やvaluation拡張を要する古典的事実
           「Eisenstein多項式は完備局所体上で全分岐拡大を与える」は
           mathlibに存在しない)——`Wₙ`の完備性(Henselian性)すら
           このセッションの仮定に含まれていないことも確認した。
      2. **原文が使う「`Ω¹_{Wₙ₊₁/Vₙ}`が`d+1`個の巡回加群の直和」という
         事実(Nakayama+elementary divisors、`hrec`の`(d+1)`という
         係数の由来)は、`pushoutKaehlerSplit`を`d+1`回反復すれば
         得られるはずと分かっているが、実際に一般の`d`について
         証明されているのは**3項(`d=2`)の場合の実証のみ**
         (`pushoutKaehlerSplit3`、line 1468近辺)——任意の`d`への
         一般化(`Fin(d+1)→Type*`で添字付けた塔の帰納法)は着手して
         いない。

      **結論(正直な評価)**: `hrec`への接続には、(1)`δn`を定義する
      2つの経路のうち**どちらか一方を独立に完成させる**(単項生成元の
      構成、または Eisenstein⟹全分岐という古典的整数論を mathlib に
      無い状態から証明する)ことと、(2)`d+1`項の`pushoutKaehlerSplit`
      の一般化、という**独立した2つの実質的な形式化タスク**が必要——
      どちらも「軽い」ものではなく、このセッション内で追加的に完成
      させるのは現実的でないと判断する。ここまでで確立した代数的な
      骨格(`falt1_theorem12_differentIdeal_length`・`falt1_
      differentIdeal_diamond_length`・`delta_tendsto_zero`・
      `pushoutKaehlerSplit3`)は、これら2つのタスクが完成すれば
      即座に組み合わさる形になっている——次回以降のセッションは
      この2つのうちどちらか一方(特に(2)の`d+1`項化は純粋に代数的で
      新しい整数論を要しないぶん、着手しやすい可能性がある)から
      始めるとよい。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **(2)の`d+1`項化に実際に着手し、想定より軽く解決できた**
      (`pushoutKaehlerSplitStep`/`pushoutKaehlerSplitBase`、commit
      `91082cdb`、`lake build`完走(169秒)・sorry無し)。`pushoutKaehler
      Split3`(3項=`d=2`のみ実証)を見て、「前の状態」を**2つの決まった
      因子ではなく任意の`Fintype`族**として抽象化すれば帰納の**1段**を
      切り出せると気づいた——`TensorProduct.piRight`(`N⊗[R](∀i,M i)
      ≃∀i,N⊗[R]M i`、任意の`Fintype`添字で無条件に成立)と
      `LinearEquiv.piCongrRight`(各成分ごとの同型をPi型全体へ持ち上げる)
      という、mathlibに**既にあった**道具を使うだけで、`pushoutKaehler
      Split3`の2段目の議論がそのまま任意個の前の因子に対して機械的に
      一般化できた。

      **設計判断**: 「`Fin(d+1)`全体を1本の定理として閉じる」ラッパー
      (`Fin.snoc`での添字の付け替えを要する)はあえて作らなかった
      ——`pushoutKaehlerSplitStep`(1個追加する帰納の1段、任意の
      `Fintype`添字の「前の状態」を受け取る)自体が、実際の使用
      (Theorem 1.2の`V_n`塔の構成、それ自体が`n`についての帰納)で
      「1段ずつ繰り返し適用する」形でそのまま使われるはずなので、
      `Fin(d+1)`版という**閉じた形での再定式化はむしろ不要**だと
      判断した——`V_n`塔を実際に構成する際にこの1段を`n`回チェイン
      すればよい。これで「(2)は純粋に代数的で新しい整数論を要しない
      ぶん着手しやすい」という前回の見立てが正しかったことが確認
      できた。

      **残る作業(更新)**: (1)の`δₙ`の実数化(単項生成元 or
      Eisenstein⟹全分岐)は依然として未解決。(2)の帰納の道具は揃った
      ので、次回は実際に`V_n`塔を`n`段構成しながら
      `pushoutKaehlerSplitStep`を適用する場面(3cの再帰的構成、
      `falt1-goal.md`の他の箇所で未着手と記録済み)に着手できる状態に
      なった。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **残る(1)も部分的に前進させた**(`falt1_theorem12_kaehler_length`・
      `falt1_theorem12_kaehler_length_eq_differentIdeal`、commit
      `69478240`、`lake build`完走(165秒)・sorry無し)。`falt1CokernelLengthEq`
      (Lemma 1.1)を`V0→Wₙ₊₁`へ直接適用してみたところ、**`Wₙ₊₁`が
      `V0`上局所的(全分岐)であることは要求しない**と判明した
      (`falt1CokernelLengthEq`自体、任意のDedekind整域`W`が`V`上
      単項生成であれば成り立つ一般論——局所性は`length_map_pow_of_
      ramificationIdx`側だけの要求だった)。つまり option (b)
      (Eisenstein⟹全分岐、mathlibに無い)は不要で、**option (a)
      (単項生成元)だけが本質的な残る前提**だと確認できた。

      `IsIntegralClosure Wₙ₊₁ V0 (FractionRing Wₙ₊₁)`も
      `falt1_theorem12_differentIdeal_length`の証明で既に使った
      `IsIntegralClosure.of_isIntegrallyClosed`で無条件に得られる
      (`Wₙ₊₁`のDedekind性・`V0`上有限であることだけで十分)。

      これで、`Wₙ₊₁`の`V0`上の単項生成元`y`(環・体どちらのレベルでも
      生成)**さえ与えられれば**、
      ```
      length_{Wₙ₊₁}(Ω¹_{Wₙ₊₁/V0}) =
        length_{Wₙ₊₁}(Wₙ₊₁/differentIdeal Wₙ Wₙ₊₁) +
        length_{Wₙ₊₁}(Wₙ₊₁/Ideal.map(algebraMap Wₙ Wₙ₊₁)(differentIdeal V0 Wₙ))
      ```
      という、Kähler微分の完全列(Exercise 13.7.4)経由と`differentIdeal`
      の塔公式経由という**2つの独立した証明経路が合流する式**が得られた
      ——これがTheorem 1.2の目標そのものに極めて近い形。

      **唯一残る条件**: `Wₙ₊₁`の`V0`上の単項生成元`y`の存在。これは
      「`Wₙ₊₁`が半局所環でも成り立つか」という問いであり、まだ
      検証していない——次回、実際の`V_n`塔構成(3c、`pushoutKaehler
      SplitStep`が使えるようになった)の中でこの生成元をどう構成する
      か、あるいは一般に存在しない場合の代替(`y`無しで済む形への
      式変形)を検討することから始める。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **上記「`y`無しで済む形」の代替を実際に見つけ、`pushoutKaehler
      SplitStep`を「長さ」レベルまで持ち上げた**(`module_length_pi_fin`・
      `module_length_pi`・`pushoutKaehlerSplitStep_length`、commit
      `1927276d`、`lake build`完走(171秒)・sorry無し)。`x`が`Wₙ₊₁`の
      `V0`上の単項生成元になれない理由(`x`の`V0`上の最小多項式の次数は
      高々`n`——`x^n = algebraMap V0 Wₙ₊₁ π`という関係が`Wₙ`を経由せず
      直接`V0`上でも成り立つため——`Wₙ₊₁`の`V0`次数`n·[Wₙ:V0]`より
      小さい)を確認し、**単項生成元を経由しない元々の原文のアプローチ
      (「`Ω¹`は`d+1`個の直和」)こそが正攻法だった**と再確認した。

      `pushoutKaehlerSplitStep`(先に確立、`LinearEquiv`)を実際に
      「長さ」の言葉に翻訳するには`Module.length`の`Pi`型への加法性
      (`Module.length R(∀i,M i)=∑i,Module.length R(M i)`)が要るが、
      mathlibには**変動する族**に対するこの一般形が見当たらなかった
      ため(`Module.length_prod`は2項、`Module.length_pi`は定数族のみ)、
      `Fin.consLinearEquiv`による`Fin n`上の帰納法+`LinearEquiv.
      piCongrLeft`(添字の付け替え)で新規に証明した。結果:
      ```
      Module.length B(Ω[B/R]) =
        (∑i,Module.length B(TensorProduct(F i)B Ω[F i/R])) +
        Module.length B(TensorProduct C B Ω[C/R])
      ```
      これで「`Ω_{Wₙ₊₁/Vₙ}`は`d+1`個の巡回加群の直和」という原文の
      主張が、**単項生成元を一切経由せず**「長さ」のレベルまで到達
      した。次回: この`pushoutKaehlerSplitStep_length`を実際の`V_n`塔
      (`V0→V1→...→Vn`、各段は単項生成——各**個別の段**は monogenic
      なので Lemma 1.1 が使える)に沿って`n`回チェインし、各因子
      `TensorProduct(F i)Wₙ₊₁ Ω[F i/Vᵢ]`(1段分の拡大)を`falt1CokernelLengthEq`
      で`differentIdeal`に変換して和を取れば、`hrec`が要求する
      `δₙ`の実数化に到達できる見込み——ここから続ける。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **`pushoutKaehlerSplitStep`をFalt1の具体的な`Wₙ₊₁`へ実際に
      適用しようとして、鍵となる既存の道具を発見した**(未完成、
      調査のみ・commit無し)。`pushoutKaehlerSplitStep`(`R B1 C B`)を
      `R:=V0,B1:=Wₙ,C:=V1(≃AdjoinRoot f),B:=Wₙ₊₁`に使うには
      `Algebra.IsPushout V0 Wₙ (AdjoinRoot f) Wₙ₊₁`(`Wₙ₊₁`が`Wₙ`と
      `V1`の`V0`上のpushout=テンソル積であること)が要る。これを
      示す道具は**このセッションより前に既に用意されていた**と判明
      した:
      - `adjoinRootTensorEquiv R C g : TensorProduct R C(AdjoinRoot g)
        ≃ₐ[C] AdjoinRoot(g.map(algebraMap R C))`(既存、`KaehlerAux.lean`
        line 1394近辺)——`R:=V0,C:=Wₙ,g:=f`とすれば`Wₙ⊗[V0]AdjoinRoot f
        ≃ₐ[Wₙ] AdjoinRoot g`(=`Wₙ₊₁`)がまさにこれ。
      - `algHomAdjoinRootOfCompat' f : AdjoinRoot f →ₐ[V0] AdjoinRoot
        (f.map(algebraMap V0 Wₙ))`(既存、line 2450近辺)——`AdjoinRoot f`
        (≃V1)から`AdjoinRoot g`(≃Wₙ₊₁)への自然な埋め込み。
      - `Algebra.IsPushout.of_equiv`(mathlib既存)——標準の
        `TensorProduct.isPushout`(`Algebra.IsPushout R C(AdjoinRoot f)
        (TensorProduct R C(AdjoinRoot f))`、無条件のinstance)を
        `adjoinRootTensorEquiv`で`AdjoinRoot g`へ運べば
        `Algebra.IsPushout V0 Wₙ(AdjoinRoot f)(AdjoinRoot g)`が出る
        ——ただし`e.toRingEquiv.toRingHom.comp(algebraMap(AdjoinRoot f)
        (TensorProduct...)) = algebraMap(AdjoinRoot f)(AdjoinRoot g)`
        という**両立性**を示す必要があり、これは`adjoinRootTensorEquiv`
        の内部定義(`Algebra.TensorProduct.tensorQuotientEquiv`・
        `tensorPolynomialAlgEquiv`・`Ideal.quotientEquivAlg`の合成)を
        `AdjoinRoot.mk g p`の形の元に対して追跡する必要がある——
        REPLで試みたが、複数層の unfold が絡み合い今回は閉じ切れな
        かった(`sorry`のまま、commitしていない)。
      **次回の具体的な着手点**: この両立性の証明(`AdjoinRoot.induction_on`
      + `Polynomial.induction_on`、または`adjoinRootTensorEquiv`
      自体に「`1⊗ₜ(AdjoinRoot.mk g p) ↦ AdjoinRoot.mk(g.map φ)(p.map φ)`」
      という`simp`補題を先に用意しておく方が近道かもしれない)を完成
      させれば、`Algebra.IsPushout V0 Wₙ(AdjoinRoot f)(AdjoinRoot g)`
      が手に入り、`pushoutKaehlerSplitStep_length`が(`ι:=Fin 1`、
      `F 0:=AdjoinRoot f`、`prev:=`自明な`Ω[Wₙ/V0]`の1点分解、という
      **最初の1段**として)Falt1の具体的な対象に対して即座に適用できる
      見込み。

      ★★2026-09-04、上記の`simp`補題(「`1⊗ₜ(AdjoinRoot.mk g p) ↦
      AdjoinRoot.mk(g.map φ)(p.map φ)`」)を実際にREPLで組み立てようと
      したところ、より詳細な技術的障害を特定した(未解決): `unfold
      ABC3.Found.Falt1.adjoinRootTensorEquiv`すると`have e1:=...;have
      hpt:=...;have hIJ:=...;e1.trans(...)`という`let`鎖の形で展開
      され、`dsimp only`でzeta簡約したあとも
      `ABC3.Found.Falt1.tensorPolynomialAlgEquiv_one_tmul`を`rw`しよう
      とすると「`↑tensorPolynomialAlgEquiv(1⊗ₜg)`(RingHomとしての
      coercion、`hIJ`内部で使われている形)」と「`tensorPolynomialAlgEquiv
      (1⊗ₜg)`(AlgEquivの直接のFunLike coercion、simp補題の形)」という
      **2種類の異なるcoercion経路**が食い違い、「target expression is
      not type-correct under the instances transparency level」という
      エラーになる。次回はこの2つのcoercionを橋渡しする補題
      (`AlgEquiv.coe_ringHom`的な、`norm_cast`用の橋渡し)を先に用意
      するか、`adjoinRootTensorEquiv`自体を(`hIJ`の証明で`tensor
      PolynomialAlgEquiv`をRingHomとして扱う代わりに)AlgEquivの
      coercionで統一する形に書き直すことから始めるとよい。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **上記の障害の真因を特定し、半分(定数の場合)を実際に解決した**
      (REPLで確認、未commit)。真因は「RingHomとしてのcoercionと
      AlgEquivの直接のcoercion」という単純な話ではなく、**`AdjoinRoot
      h`と`Polynomial _ ⧸ Ideal.span{h}`が定義的には同じ型でも、
      別々に定義された`Semiring`/`Algebra`インスタンス
      (`AdjoinRoot.instAlgebra`系 vs `Ideal.Quotient.ring`/
      `Ideal.instAlgebraQuotient`系)を持つ、という mathlib 自体の
      構造的な instance diamond**だった(`unfold`すると`Ideal.quotientEquivAlg`
      が`Ideal.Quotient`側のインスタンスで型付けされた項を返すが、
      ゴールは`AdjoinRoot`側のインスタンスを要求する、という
      「Application type mismatch」で確認)。

      **回避策**: `unfold`して内部を直接いじるのではなく、
      `adjoinRootTensorEquiv g`を**ブラックボックスのAlgEquivとして
      扱い**、高レベルAPI(`AlgEquiv.commutes`・`.restrictScalars`)
      だけで両立性を示す。`Polynomial.ringHom_ext'`(2つの`Polynomial
      R →+* S`が`C`(定数)と`X`(生成元)の像で一致すれば等しい)で
      問題を「定数の場合」と「`X`の場合」に分割し、**定数の場合は
      完全に解決した**:
      ```
      adjoinRootTensorEquiv g (1⊗ₜ algebraMap R(AdjoinRoot g) r)
        = algebraMap R(AdjoinRoot(g.map φ)) r
      ```
      は`AlgEquiv.commutes((adjoinRootTensorEquiv g).restrictScalars R) r`
      (`unfold`不要、instance diamond を回避)+`1⊗ₜ(r•1) = algebraMap
      R(TensorProduct...) r`という小さな計算(`TensorProduct.tmul_smul`・
      `Algebra.TensorProduct.one_def`)だけで閉じた。

      **残る`X`の場合**(生成元自体の像、`adjoinRootTensorEquiv g
      (1⊗ₜ root g) = root(g.map φ)`)は、`AlgEquiv.commutes`のような
      「algebraMapの像なら自動的に分かる」経路が使えない(`root`は
      algebraMapの像ではない)ため、依然として`unfold`(=instance
      diamond)を避けられず未解決。次回はこの`X`の場合を、
      (a) instance diamond を明示的に`convert`/`cast`で橋渡しする、
      (b) `AdjoinRoot`と`Ideal.Quotient`のインスタンスが実際に等しい
      という一般的な橋渡し補題をまず用意する、のどちらかで攻める
      ことから始める。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04、
      **(a)(b)どちらでもない第3の道で、点ごとの挙動の証明を完全に
      解決した**(`adjoinRootTensorEquivFwd`・`adjoinRootTensorEquivFwd_
      one_tmul_mk`・`eval2_root_eq_mk`、commit `6740e740`、`lake build`
      完走(171秒)・sorry無し)。`adjoinRootTensorEquiv`自体を修正
      しようとするのをやめ、**`Ideal.Quotient`のAPIを一切経由しない
      別経路の線形写像**を`AdjoinRoot`ネイティブAPIだけで新規に構成
      した:
      - `adjoinRootTensorEquivFwd g := TensorProduct.AlgebraTensorModule.
        lift`で`(c,x) ↦ c • algHomAdjoinRootOfCompat' g x`という双線形
        写像を持ち上げただけ(`algHomAdjoinRootOfCompat'`は`AdjoinRoot.map`
        ベースの既存の道具、instance diamond を経由しない)。
      - その点ごとの挙動(`1⊗ₜ(AdjoinRoot.mk g p) ↦ AdjoinRoot.mk(g.map φ)
        (p.map φ)`)は、`AdjoinRoot.lift_mk`(`AdjoinRoot.map`の内部構成)
        + 新規補題`eval2_root_eq_mk`(`eval₂(AdjoinRoot.of q)(root q)h =
        AdjoinRoot.mk q h`、`Polynomial.induction_on`による3行の帰納法)
        だけで**完全に**証明できた——`X`の場合(生成元自体)を含めて
        すべて解決している(`Polynomial.induction_on`の`C`・`add`・
        `monomial`の3ケースがそのまま定数・和・生成元の場合をカバー
        するため、`ringHom_ext'`のような分割すら不要だった)。

      **これで`adjoinRootTensorEquiv`の代わりに`adjoinRootTensorEquivFwd`
      を使えば、instance diamond を一切経由せず`Algebra.IsPushout`へ
      接続できる見込みが立った。** 残るのは`adjoinRootTensorEquivFwd`
      の**全単射性**のみ(まだ未証明)——`g`がmonicなら両辺は`C`上
      階数`deg g`の自由加群になるはずなので、(全射性を`root(g.map φ)`
      の冪が像に入ることから示し)+(有限自由加群間の全射線形写像は
      階数が等しければ自動的に単射、という一般論)で閉じられる見込み。
      次回はここから続ける。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続)、
      **`adjoinRootTensorEquivFwd`の全単射性が完成した**(階数比較の
      迂回路ではなく、**明示的な逆写像を直接構成する**という、もっと
      素朴だが完全に閉じる経路で解決):
      - `adjoinRootTensorEquivInv g := AdjoinRoot.lift (algebraMap C
        (R⊗C(AdjoinRoot g))) (1⊗ₜroot g) _`——`root(g.map φ) ↦ 1⊗ₜroot g`・
        `C ↦ algebraMap C _`となる`AdjoinRoot.lift`。well-definedness
        (`g.map φ`が`1⊗ₜroot g`で消える)は`Polynomial.aeval_algHom_apply`
        (`aeval`のAlgHom越しの自然性、`1⊗ₜ-`を`includeRight`として見る)
        +`AdjoinRoot.aeval_eq`(`aeval(root f)f=mk f f=0`)だけで閉じた
        (ここでも`Ideal.Quotient`は一切経由しない)。
      - 点ごとの挙動を2つの補題(`adjoinRootTensorEquivInv_monomial`
        =`C a*X^n`での挙動、`adjoinRootTensorEquivInv_mk_map`=`mk(g.map φ)
        (p.map φ)`での挙動)で確定し、`C`線形性を`adjoinRootTensorEquivInv_
        smul`(`algebraMap`との両立、`AdjoinRoot.lift_of`)で別途確立。
      - **往復1**(`fwd∘inv=id`、`adjoinRootTensorEquiv_roundtrip1`):
        `AdjoinRoot.induction_on`→`Polynomial.induction_on`(`C`・`add`・
        `monomial`)の二重帰納法。`monomial`のケースは`pow_succ`で`X^(n+1)
        =X^n*X`に分解し、`adjoinRootTensorEquivInv_monomial`(`n=1`の
        特殊形)と`adjoinRootTensorEquivFwd_one_tmul_mk`を`Algebra.
        TensorProduct.tmul_mul_tmul`で繋いだ。
      - **往復2**(`inv∘fwd=id`、`adjoinRootTensorEquiv_roundtrip2`):
        `TensorProduct.induction_on`→`AdjoinRoot.induction_on`の二重
        帰納法。`tmul`ケースは`c⊗ₜz = c•(1⊗ₜz)`という`TensorProduct.
        smul_tmul'`による書き換えで`C`線形性(`adjoinRootTensorEquivInv_
        smul`)に帰着させ、`1⊗ₜ-`の場合だけ`adjoinRootTensorEquivFwd_
        one_tmul_mk`+`adjoinRootTensorEquivInv_mk_map`の合成で閉じた。
      - **`adjoinRootTensorEquivFwd_bijective`**: 上記2つの往復から
        `Function.bijective_iff_has_inverse`で直ちに従う。

      これで`AdjoinRoot`/`Ideal.Quotient`のinstance diamondに一度も
      触れずに、`TensorProduct R C (AdjoinRoot g) ≃ₗ[C] AdjoinRoot(g.map φ)`
      という**完全な線形同型**(の存在、`LinearEquiv.ofBijective`経由)
      が手に入った。次は`Algebra.IsPushout.of_equiv`(または`IsBaseChange.
      of_equiv`)に`adjoinRootTensorEquivFwd_one_tmul_mk`から従う
      tmul整合性条件を渡し、`Algebra.IsPushout V0 Wn (AdjoinRoot f)
      (AdjoinRoot g)`(=具体的な`Falt1`のオブジェクトでのpushout性)を
      確立、`pushoutKaehlerSplitStep`/`_length`を実際の`V_n`/`W_n`の塔に
      接続する。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続々)、
      **`Algebra.IsPushout`への接続そのものが完成した**(commit分は
      次項で記録)。決定的な近道が判明: `Algebra.IsPushout R S R' S'`は
      mathlibで単に`IsBaseChange S (IsScalarTower.toAlgHom R R' S').
      toLinearMap`のラッパー(唯一のフィールド`Algebra.IsPushout.out`)
      だと`#print`で確認できた——つまり`TensorProduct.isPushout'`
      経由で「既知のpushout」を`of_equiv`で移送する回り道は不要で、
      `adjoinRootTensorEquivFwd_bijective`から**直接**`IsBaseChange`を
      経て`constructor`一発で`Algebra.IsPushout`に落ちる:
      - `falt1_isBaseChange_adjoinRoot`: `IsBaseChange.of_equiv`に
        `LinearEquiv.ofBijective adjoinRootTensorEquivFwd
        adjoinRootTensorEquivFwd_bijective`を渡し、tmul整合性条件を
        `adjoinRootTensorEquivFwd_one_tmul_mk`+新規補題
        `algHomAdjoinRootOfCompat'_mk`(`algHomAdjoinRootOfCompat'`の
        `AdjoinRoot.mk`での挙動)で確認。
      - `falt1AdjoinRootAlgebra`: `algHomAdjoinRootOfCompat' g`
        (=`V1→W1`の橋渡し写像)を`RingHom.toAlgebra`で
        `Algebra(AdjoinRoot g)(AdjoinRoot(g.map φ))`インスタンスに
        変換。`instance`ではなく`def`(`@[reducible]`)のまま保った——
        `Falt1`の具体的オブジェクトへ適用する際、既存の`Algebra V1 Wn1`
        インスタンスとの衝突(前回セッションで診断済みの問題)を
        避けるため。
      - `falt1AdjoinRoot_isScalarTower`: `IsScalarTower.of_algebraMap_eq`
        +`algHomAdjoinRootOfCompat'`の`.commutes`(`R`上の`AlgHom`である
        こと)から直ちに従う。
      - `falt1_isPushout_adjoinRoot`(到達点): `constructor`+`show`+
        `falt1_isBaseChange_adjoinRoot`の3行で`Algebra.IsPushout R C
        (AdjoinRoot g)(AdjoinRoot(g.map φ))`が確立。`AdjoinRoot`/
        `Ideal.Quotient`のinstance diamondには最後まで一度も触れて
        いない。

      これで`pushoutKaehlerSplitStep`/`_length`が要求する`[Algebra.
      IsPushout R B1 C B]`インスタンスに、`B1=AdjoinRoot g`(=`V1`)・
      `B=AdjoinRoot(g.map φ)`(=`Wn1`)として直接対応する一般形が
      手に入った。残る仕事は、この抽象形を具体的な`Falt1`の`V_n`・
      `W_n`の塔(`falt1AdjoinRootGeneratorFull`等で既に構成済みの
      オブジェクト)に**実際にインスタンス化**すること——ここで
      「`Algebra V1 Wn1`インスタンス衝突」(前回セッションで診断済み)
      に再度注意しながら進める必要がある。

      ★2026-09-04(続々々)、**インスタンス化を試みる前に、2つの構造的な
      ギャップが判明した**(`falt1_cancelConductorDelta_assembled`の
      シグネチャを再確認して発見):
      1. **base ringのずれ**: `falt1_isPushout_adjoinRoot`は`g :
         Polynomial R`(`R`は素の底環)についての命題だが、実際の
         `V1 := integralClosure V0 (AdjoinRoot fK)`・`Wn1 :=
         integralClosure Wn (AdjoinRoot gK)`は、`fK`・`gK`が
         `FractionRing V0`・`FractionRing Wn`(**分数体**)上の多項式
         として定義されている(Eisenstein性の議論に分数体が必要な
         ため)。よって`falt1_isPushout_adjoinRoot`を使うには
         `R = FractionRing V0`・`C = FractionRing Wn`とせざるを得ず、
         得られるのは`Algebra.IsPushout (FractionRing V0)
         (FractionRing Wn) (AdjoinRoot fK) (AdjoinRoot gK)`であって、
         `pushoutKaehlerSplitStep`が要求する`Algebra.IsPushout V0 Wn
         V1 Wn1`(=**整数環**のレベルでのpushout)ではない。
      2. **`AdjoinRoot`と`integralClosure`のずれ**: `V1`・`Wn1`は
         `AdjoinRoot fK`・`AdjoinRoot gK`**そのもの**ではなく、その中の
         `V0`・`Wn`の**整閉包**(=分岐拡大の「整数環」)である——
         `AdjoinRoot fK`は`FractionRing V0`上有限次数の**体**の
         拡大体そのものであり、`V1`はその中の真部分環。
      3. **数学的な疑義(最重要)**: 2026-09-04(続)に記録した
         Brinon-Conrad Exercise 13.7.4 の発見(「`Ω¹_{B1/A}`の生成元の
         個数は second fundamental exact sequence + Nakayama で出す」
         という手法)と照らすと、Faltings自身の証明は**`Wₙ₊₁`が`Vₙ₊₁`
         (=`V1`)と`Wₙ`の pushout である、という同型そのものは
         主張していない**可能性が高い——`Wₙ₊₁`は非正規かもしれない拡大の
         整数環であり、pushout(=テンソル積)は一般には整閉ではない
         ため、`Algebra.IsPushout V0 Wn V1 Wn1`が字義通り成り立つとは
         限らない。むしろ Faltings は「`Ω¹_{Wₙ₊₁/Vₙ₊₁}`を`d+1`個の
         生成元をもつ加群として直接扱う」形で pushout の**厳密な同型**
         を回避している可能性がある。

      ★結論: `falt1_isPushout_adjoinRoot`(分数体レベルでの一般公式)は
      それ自体正しく完成した独立した成果だが、これを`V1`/`Wn1`
      (整数環レベル)へ**そのまま**インスタンス化するのは早計——
      次回はまず「`Algebra.IsPushout V0 Wn V1 Wn1`は本当に成り立つか」
      を(反例や、成り立つとしても追加の仮定(例えば`Wn`が`V0`上
      不分岐、など)が要るかを含めて)再検討することから始める。
      性急に誤った instantiation を試みるより、この構造的な疑問を
      先に解消する方が着実。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **上記の疑義を、Theorem 1.2 の証明原文(物理p.5、上記の逐語訳)を
      再読して確定させた——`Algebra.IsPushout V0 Wn V1 Wn1` は
      文字通り成り立たない、が確定。むしろ本来の適用先が判明した。**

      原文第2文からそのまま読める:「`Ω_{W_n/V_n}⊗W_{n+1} →
      Ω_{W_{n+1}/V_n}`」という**写像**を考え、その kernel・cokernel を
      `p^{δ_n-δ_{n+1}}` で評価する、というのが Theorem 1.2 の証明の
      **核心そのもの**——`W_n⊗_{V_n}V_{n+1}`(=もし pushout が成り立てば
      `W_{n+1}`と同じになるはずの対象)と実際の`W_{n+1}`(=正規化後の
      整数環)との**ズレ**(`p^{δ_n-δ_{n+1}}`で零化される)を測ることが
      証明の仕事であり、このズレが**ちょうど0**(=`Algebra.IsPushout`
      成立)なら証明の大部分は不要になってしまう。ゆえに
      `Algebra.IsPushout V0 Wn V1 Wn1` は Theorem 1.2 の**主張が正しい
      限り一般には偽**(反例探しは不要——原文自身がそれを否定する
      形で構成されている)。

      **一方、`pushoutKaehlerSplitStep`/`falt1_isPushout_adjoinRoot`の
      本来の適用先は`W_n→W_{n+1}`の橋ではなく、`V_n→V_{n+1}`という
      "同時添加"の側だった**(3a/3b の作業ログに既に記録済みの方針を
      再確認・確定):「典型例」段落(`V_{n+1}` は `T_i` の `p` 乗根
      `d` 個 + `1` の `p^{n+1}` 乗根 `1` 個の**同時添加**)が示す通り、
      `V_{n+1}` 自体が(`V_n`上の)`d+1` 個の単一生成拡大の
      **反復pushout**として——`V_{n+1} = AdjoinRoot g_1 ⊗_{V_n} ⋯
      ⊗_{V_n} AdjoinRoot g_{d+1}`という形で——**定義そのものによって**
      構成できる(ここには正規化のズレが原理的に無い、pushout性は
      "証明すべき事実"ではなく"定義から従う自明な事実"になる)。
      `falt1_isPushout_adjoinRoot`は`R`が任意の`CommRing`で
      成り立つ(`FractionRing`を要求しない、`g:Polynomial R`のまま)
      ため、`R=V_n`・各段の`C=AdjoinRoot g_i`として**まさにこの
      構成に直接使える**——これが今回発見した本来の使い道である。

      **次に必要な新規構築**(未着手、次回以降の課題として明記):
      1. `V_{n+1}`を`d+1`個の単一生成`AdjoinRoot`因子の反復
         テンソル積として具体的に定義する(`pushoutKaehlerSplitStep`
         を`d+1`回適用する帰納、3aで既に汎用化済みの1段をそのまま
         使う)。
      2. この`V_{n+1}`が実際に(Faltings の意図通り)整閉であること
         ——各因子が`V_n`上不分岐または適切な独立性を持つ場合に
         テンソル積が整閉であるという古典事実(3b で既に記録した
         「conductor-discriminant」の周辺、`Algebra.discr`系)を要する
         可能性が高い。
      3. `Ω_{V_{n+1}/V_n}⊗_{V_{n+1}}W_{n+1}`(原文が実際に使う対象、
         `V_{n+1}`側の分解を`W_{n+1}`まで底変換したもの)への接続——
         `pushoutKaehlerSplitStep_length`が`Ω_{V_{n+1}/V_n}`自体の
         分解を与えるので、そこから`⊗_{V_{n+1}}W_{n+1}`する底変換は
         別途(`TensorProduct.AlgebraTensorModule`系)。

      この経路は「`W_n`/`W_{n+1}`のinstance衝突を避けながら
      `Algebra.IsPushout V0 Wn V1 Wn1`を確立する」という(誤った)
      目標を追うより筋が良い——`V_n`側だけで完結する構成であり、
      `Wn`は最後の底変換の段階まで一切登場しない。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **上記1(反復pushoutで`V_{n+1}`を構成)の**核となる2つの道具**を
      完成させ、1段の反復pushoutが実際に繋がることを確認した**
      (commit `4943ba71`):
      - `adjoinRootTensorAlgHom`/`_eq_fwd`/`_bijective`/
        `adjoinRootTensorAlgEquiv`: `adjoinRootTensorEquivFwd`(線形
        同型)を`AlgEquiv`として再構成(`Algebra.TensorProduct.lift`
        で`AlgHom`化し、pure tensorでの一致から全単射性を移送)——
        `mapBaseChange_injective_transport`のような`AlgEquiv`を要求
        する道具に、instance diamond を再導入せず渡せるようになった。
      - `mapBaseChange_injective_adjoinRoot_direct`:
        `mapBaseChange_injective_adjoinRoot_tensor`(`TensorProduct`版)
        を上の`AlgEquiv`で`AdjoinRoot(g.map φ)`側へ移送し、
        `pushoutKaehlerSplitStep`の`hinj`を Lemma 1.1 型の非零因子
        条件1つに直接帰着させた。

      実際に`R`(=`V_n`役)・`g1 g2 : Polynomial R`を使い、`B1:=R`
      (`pushoutKaehlerSplitBase`)から`falt1_isPushout_adjoinRoot g1`+
      `mapBaseChange_injective_adjoinRoot_direct g1`で`pushoutKaehler
      SplitStep`を1回適用できることをREPLで確認した(`letI`で
      `falt1AdjoinRootAlgebra`/`falt1AdjoinRoot_isScalarTower`を
      instanceとして局所登録するだけで、残りの`Algebra R B`・
      `IsScalarTower`系はすべて`infer_instance`で自動解決した)。

      **2段目(`g2`をさらに添加)を試みたところ、別の——3aで既に
      予告されていた——障害に到達した**: 1段目の出力
      `(∀i:Fin1,F i) × TensorProduct(AdjoinRoot g1)B2 Ω[AdjoinRoot
      g1/R]`を2段目の`prev`が要求する`∀i:ι',F' i`という**1つの
      Pi型**に組み替える必要があるが、これを`Option ι`(または
      `Fin(n+1)`)で束ねようとすると、`Ω[(Option.elim i C F)/R]`と
      いう**型`i`に依存する式そのもの**が`CommRing(Option.elim i C
      F)`のようなinstanceを要求し、これは`i`が自由変数である限り
      (`i=none`/`some j`と具体的に分岐できないため)`infer_instance`
      では解決できないことを実測で確認した(3aのdocstringが予告して
      いた「`recursive Type-valued def with dependent instances`」
      という地雷に、より具体的な形で到達した)。

      **今後の対処方針**(次回、着手する場合の設計指針として記録):
      `Option.elim`/`Fin.snoc`で型族を都度組み替えるのではなく、
      3aのdocstringが示唆した通り**Σ束ね**(`structure RAlg (R)
      where carrier : Type*; [ring : CommRing carrier]; [alg :
      Algebra R carrier]`)を使い、「`i`番目の環+その環自身が運ぶ
      instance」を**1つのデータ**として扱う設計にすれば、
      `(Option.elim i c f).ring`のような**射影**は`i`の分岐無しに
      常に型付けされる(Σ型の外側の`Option.elim`の結果自体は自由変数
      `i`のままでも、`.ring`フィールドの射影は無条件に使える)ため、
      この instance 解決の壁を回避できる見込み。ただし`Ω[·/R]`
      (`KaehlerDifferential`)をΣ束ねの`carrier`越しに使うための
      追加の配線(`Algebra R (Option.elim i c f).carrier`をΣの
      フィールドとして運ぶ、等)が必要で、まだ試していない。

      ★2026-09-05、**上記のΣ束ね方針を実際にREPLで検証し、核となる
      部分は機能することを確認した**:
      ```lean
      structure RAlg (R : Type*) [CommRing R] where
        carrier : Type*
        [ring : CommRing carrier]
        [alg : Algebra R carrier]
      attribute [instance] RAlg.ring RAlg.alg
      ```
      とすると、`(C F : RAlg R)`・自由変数`i : Option ι`に対して
      `CommRing (Option.elim i C F).carrier`・`Algebra R (Option.elim
      i C F).carrier`が**`infer_instance`で即座に解決**し、
      `Ω[(Option.elim i C F).carrier⁄R]`という型そのものも問題なく
      elaborate できることを確認した——予想通り、射影(`.ring`/`.alg`
      フィールド)は`i`の分岐無しに常に使えるため、地雷を回避できる。

      ただし`pushoutKaehlerSplitStep`が要求する残りのinstance
      (`∀i,Algebra(F i)B1`・`∀i,IsScalarTower R(F i)B1`等、**現在の
      累積環`B1`/`B`との関係**を述べるもの)は、`RAlg`が`R`上の
      代数であることしか記録していないため**この束ねだけではカバー
      できない**——`B1`/`B`は帰納の段ごとに変わる対象であり、
      「`i`番目の環が現在の累積環にどう埋め込まれるか」という情報は
      別途(各段で`Fin.cases`的に具体的な`i`へ場合分けして手で
      instance登録する、または`RAlg`をさらに拡張して「今までの塔への
      整合的な埋め込み系」まで一緒に運ぶ、のいずれか)必要になる——
      後者は事実上「圏論的な帰納系(directed system)」の随伴物になり、
      本腰を入れた設計判断が要る規模の作業と判断し、ここでは着手を
      保留する(次回以降の課題として明記のみ)。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **「後者」(`RAlg`を現在の累積環`B`への埋め込みごと運ぶ拡張)を
      実際に実装し、`pushoutKaehlerSplitStep`が要求する`B1`関係・`B`
      関係のinstance両方を自由変数の添字`i`のままで供給できることを
      確認した**(`RAlg`/`RAlgOver`/`RAlgOver.lift`/`RAlgOver.lift_
      isScalarTower`の4点、`lake build`完走・sorry無しで
      `Found/Falt1/KaehlerAux.lean`に追加済み):
      - `RAlgOver (R T)`: `RAlg R`に加え、具体的な「現在の累積環`T`」
        への埋め込み(`Algebra carrier T`・`IsScalarTower R carrier
        T`)まで一緒に運ぶΣ束ね。
      - `RAlgOver.lift`: `B1`への埋め込みを、`B1↪B`という1段の拡大
        越しに`B`への埋め込みへ運ぶ(合成`carrier→B1→B`を
        `RingHom.toAlgebra`で新しい`Algebra carrier B`として登録)。
      - `RAlgOver.lift_isScalarTower`: 運んだ埋め込みが
        `carrier→B1→B`という`IsScalarTower`をなすこと。

      **鍵となった実装上の教訓**(次回、この種の「Σ束ねの値をtermで
      運ぶ」パターンを書く際に直接使える): `letI`(値の透明性を保つ)
      と`haveI`(値を消してしまう、命題には無害だが**データ**の
      instanceには不可)の違いが決定的だった——`haveI hAlgB := fun i
      => (F i).lift.algT`のように`haveI`で登録すると、後で`(F i).
      lift.towerT`のような**別の場所から来る同じデータ**と`hAlgB i`
      が(命題上は等しくても)**構文的に区別**されて`Type mismatch`に
      なる。`letI`に変えるだけで、透明な`let`束縛として`hAlgB i`が
      `(F i).lift.algT`へ展開可能になり、この不一致がすべて消えた。

      **残る未完成部分**(次回への持ち越し、正確に特定できている):
      上記4点を実際に`Option ι`のPi型全体へ組み立てる最後の1段
      (`LinearEquiv.piOptionEquivProd`との合成、`pushoutKaehlerSplit
      StepOption`という名を予定)で、`none.elim`/`(some i).elim`が
      `rfl`で潰れるはずの`Option.elim i C F`絡みの式が`Type mismatch`
      (instanceスロットの不一致)を起こす。原因はおそらく
      `LinearEquiv.piOptionEquivProd`の`M`が**暗黙引数**であることに
      気付かず最初は位置引数で渡していた(`(M := ...)`と明示して
      解決した)のと同種の、まだ特定しきれていない別の暗黙引数の
      取り違えである可能性が高い——次回は`set_option pp.explicit
      true`等でinstance引数を完全に表示させてから`convert`で1つずつ
      潰す、という方針で再開すること。

      ★2026-09-05、**上記の診断を1段深く進めた**——原因は
      `set_option pp.explicit`で見える表面的な引数の取り違えではない
      ことを確認した。切り分けの結果:
      - `(Option.elim none C F).carrier = C.carrier`(`none`のケース)・
        `∀i,(Option.elim(some i)C F).carrier=(F i).carrier`(`some`の
        ケース)は、`RAlg`/`RAlgOver`でラップした**同じ設定**で
        **単独では**すべて`rfl`で通る(積の型・Pi型としてまとめても
        `rfl`で通る)ことを確認した。
      - ところが`LinearEquiv.piOptionEquivProd B (M:=...)).symm`を
        `have h := ...`で受けた後、`exact h`・`cast (by rfl) h`・
        `dsimp only [] at h`(空simpによる強制簡約、"no progress"の
        結果)のいずれも`h`の型と目標の型を一致させられない——つまり
        **単独の成分はすべてrflで潰れるのに、`LinearEquiv`として
        束ねた瞬間に一致しなくなる**。
      - これは「表面的な暗黙引数の取り違え」ではなく、
        `LinearEquiv.piOptionEquivProd`自身の内部実装が(`M none`・
        `M(some i)`を経由する際に)**目に見えない場所で追加の
        instance変換**(例えば`Module`インスタンスの正規化や
        `AddCommMonoid`側の別経路)を挟んでいる可能性を示唆する——
        `set_option pp.all true`での完全ダンプは1780行に達し、
        このセッション内での特定は断念した。
      - 次回の再開方針(更新): `pp.explicit`ではなく、
        `LinearEquiv.piOptionEquivProd`の**ソース定義**
        (mathlib`Mathlib/Algebra/Module/Pi.lean`または近傍)を直接
        読み、`M none`/`M(some i)`がどう扱われているかを確認する
        ところから始める方が近道の見込み。あるいは`piOptionEquivProd`
        を使わず、`Equiv.piOptionEquivProd`(型レベル)+
        `LinearEquiv.ofLinear`で自前に組み立て直す方が、
        instance変換の余地が無く安全かもしれない。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **原因を`LinearEquiv.piOptionEquivProd`のソース
      (`Mathlib/LinearAlgebra/Pi.lean:527`)まで実際に読み、当たりを
      さらに絞り込んだ**:
      ```lean
      def piOptionEquivProd {ι} {M : Option ι → Type*} [...] :
          ((i : Option ι) → M i) ≃ₗ[R] M none × ((i : ι) → M (some i)) :=
        { Equiv.piOptionEquivProd with map_add' := ...; map_smul' := ... }
      ```
      で、`Equiv.piOptionEquivProd`(`#print`で確認)の`invFun`は
      `fun x a => Option.casesOn a x.1 x.2`——**`Option.elim`ではなく
      `Option.casesOn`**を使っている。そこで`Option.casesOn`に
      統一して(`motive`を明示して)再試行したが、**それでも同じ
      `Type mismatch`が残った**——つまり`Option.elim`対`Option.
      casesOn`という表面的な構文の違いが原因ではなかった。

      さらに、`Fintype ι`・`DecidableEq ι`を含めて元の状況を完全に
      再現した上でも「`(Option.casesOn none C F).carrier = C.carrier`」
      という**型の等式そのもの**は単独では`rfl`で通ることを再確認
      した。これは、**型レベルでは一致しているのに、`LinearEquiv`
      として束ねると一致しなくなる**ことを意味する——`≃ₗ[B]`は
      `Module B`インスタンスを暗黙に運ぶため、`TensorProduct C.carrier
      B Ω[C.carrier/R]`の`AddCommMonoid`/`Module B`インスタンスが、
      (a)私の文脈で直接得られる経路と(b)`piOptionEquivProd`内部で
      `M none`として間接的に得られる経路とで、**型としては等しくても
      instanceとしては異なる経路で解決されている**(instance
      diamond)可能性が高い、という診断に至った。

      **これは今回のセッション序盤に発見・解決した`AdjoinRoot`/
      `Ideal.Quotient`のinstance diamond(`adjoinRootTensorEquivFwd`
      節参照)と同種の構造の問題**であり、あのときと同じ処方箋
      (「diamondと戦うのではなく、自前で新しい写像を組み立てて
      diamondを経由しない」)が有効な見込みが高い——次回は
      `piOptionEquivProd`を再利用しようとするのをやめ、`Equiv.
      piOptionEquivProd`(型レベルの土台)の上に**自分の文脈の
      instanceだけ**を使って`map_add'`/`map_smul'`を直接証明する、
      という`adjoinRootTensorEquivFwd`と同じ戦略で組み立て直す
      ことから始めるとよい。

      ★実際にこの「自前で組み立て直す」経路を試みたところ(`refine
      LinearEquiv.mk (toLinearMap := LinearMap.mk (AddHom.mk ...) ...)
      ...`を1つの式で書く形)、依存的な`match`を含む式を1つの`refine`
      に詰め込みすぎたことが原因と見られるuniverse制約の停留
      (`stuck at solving universe constraint`)に遭遇した——
      `adjoinRootTensorEquivFwd`のときのように**小さな`have`/`set`を
      1つずつ積む**(1つの巨大な`refine`にしない)形へ分解すれば
      解決できる見込みが高いが、このセッションでは完了できなかった。
      次回はここ(`prodOptionPiEquivFresh`という仮の名で試作開始
      済み)から、`toFun`・`map_add'`・`map_smul'`・`invFun`・
      `left_inv`・`right_inv`を**別々の`have`として先に確立してから**
      `LinearEquiv.mk`に渡す、という手順で再開すること。

      ★実際にこの手順(`have`を個別に確立→`LinearEquiv.mk
      (LinearMap.mk (AddHom.mk toF hadd) hsmul) invF hleft hright`)
      で組み立て直したところ、**anonymous constructorの`Prod`取り違え
      は解消した**(`LinearEquiv.mk`という名前付きコンストラクタを
      使うことで回避できた——`{...}`記法固有の問題だったと判明)。
      しかし残る`Type mismatch`は解消しなかった——しかも今回は
      **両辺の`pp`出力が文字どおり完全に同一**(1文字違わず一致)な
      のに`Type mismatch`になる、という決定的な事実を確認した。
      これは**pretty printerには表示されないレベルでの instance の
      不一致**(`Module B (∀i,M i)`のような合成インスタンスが、
      「私が構成した項」と「目標の型」とで異なる合成経路を通って
      いる)であることの動かぬ証拠——`Option.elim`/`Option.casesOn`
      の構文的な違いという当初の仮説は完全に排除された。次回は
      `set_option pp.all`の出力(1780行、今回は読み切れず)を
      `offset`/`limit`で分割して丹念に比較するか、あるいは
      `Module B (∀i,M i)`の2つの経路をそれぞれ`(inferInstance :
      Module B _)`で個別に取り出し`Eq`ではなく`HEq`や`Subsingleton`
      経由で橋渡しする、という方針から始めること。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **真の原因が判明し、`pushoutKaehlerSplitStepOption`が完全に
      完成した**(commit分は次項参照)。instance diamondではなかった
      ——原因は2つ、どちらも今回のセッションで初めて可視化できた:

      1. **`RAlg`/`RAlgOver`の`carrier : Type*`が明示的なuniverse
         変数を持たなかったこと**。`Type*`は出現ごとに独立した
         universeメタ変数を生成するため、`F : ι → RAlgOver R B`と
         別の場所での`RAlgOver R B`の使用とで、`carrier`のuniverseが
         食い違いうる。`(inferInstance : Module B (∀i, ...))`を
         `Pi.module`経由の項と比較する診断コードを書いたところ、
         「`F`の型`ι → RAlgOver.{u_1,u_2,u_5}R B`が期待される型
         `ι → RAlgOver.{u_1,u_2,u_4}R B`と一致しない」という
         **universe添字の食い違いそのもの**がエラーとして出た
         (`set_option pp.universes true`で可視化)。修正:
         `universe uFalt1R uFalt1Car uFalt1T uFalt1B`を明示的に
         宣言し、`RAlg (R : Type uFalt1R)`・`RAlgOver (R : Type
         uFalt1R)(T : Type uFalt1T)`・`carrier : Type uFalt1Car`と、
         使用箇所すべてで`RAlgOver.{uFalt1R, uFalt1Car, uFalt1T}`の
         ように**明示的にuniverse引数を揃える**。
      2. **さらに深い、真の最終ブロッカー**: 上記のuniverse修正を
         入れても`Type mismatch`は消えなかった——理由は
         `(A × B) ≃ₗ[B] C`という式を書く際、**外側の括弧を省略する
         と`A × (B ≃ₗ[B] C)`に構文解析される**(`×`が`≃ₗ[·]`より
         結合が弱い/優先順位が異なる)という**単純な構文の罠**
         だった。実際、`pp.universes`で完全一致するように直した後の
         `Type mismatch`のエラーメッセージをよく読むと、「私の項の
         型」は`(A×B)≃ₗC`なのに「期待される型」は`A×(B≃ₗC)`(`Prod`
         が一番外側)になっていた——これがセッション全体を通して
         「表面上は同一なのに`Type mismatch`」に見えていた**真の
         原因**。修正は単純: `((TensorProduct C.carrier B Ω[C.carrier
         ⁄R]) × (∀i,...))≃ₗ[B](...)`と**外側の積全体を明示的に
         括弧で囲む**だけ。

      これで`prodOptionPiEquivFresh`(のちに`pushoutKaehlerSplit
      StepOption`に統合)が`have`/`set`の個別確立
      (`toF`・`hadd`・`hsmul`・`invF`・`hleft`・`hright`)+
      `LinearEquiv.mk (LinearMap.mk (AddHom.mk toF hadd) hsmul) invF
      hleft hright`(`{...}`匿名コンストラクタは`Prod`側と誤って
      マッチする別のバグがあったため、名前付きコンストラクタで回避)
      で完全に`lake build`を通った。`pushoutKaehlerSplitStepOption`
      は`pushoutKaehlerSplitStep`の出力(積の形)を次段の`prev`に
      そのまま使える`Option ι`のPi型へ変換する——これで`d+1`個の
      因子を`n`回繰り返し適用してVₙ₊₁を構成する、という3aで構想した
      道が**工学的には完全に開通した**。

      ★続けて`pushoutKaehlerSplitStepOption_length`(長さ版)も
      即座に完成させた(commit分は次項)——`module_length_pi`
      (`Option ι`添字のPi型にも無条件成立)+`Fintype.sum_option`
      (`∑i:Option ι,f i = f none + ∑j:ι,f(some j)`)を組み合わせる
      だけで、`Module.length B(Ω[B/R]) = Module.length B(新因子の
      寄与) + ∑ⱼ Module.length B(前段の各因子の寄与)`という、
      `pushoutKaehlerSplitStep_length`の`Option ι`連鎖版が得られた。
      これで反復pushoutを`n`回適用した結果の**長さ**も、各段で
      追加した因子の長さの和として追跡できる基盤が揃った。

      ★さらに、不要な仮定`[Algebra C.carrier B1]`(`pushoutKaehler
      SplitStep`自身は要求しておらず、証明本文でも一度も使っていな
      かった——新しい因子`C`が前段の累積環`B1`に埋め込まれる必要は
      無い、むしろ一般には逆方向)を発見・削除した(commit分は次項)。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **`AdjoinRoot`因子2個の反復pushoutを実際にend-to-endで構築した**
      (`falt1_pushoutKaehlerSplitStepOption_adjoinRoot_example`、
      `pushoutKaehlerSplit3`と同じ役割の小さな具体例)。`falt1_
      isPushout_adjoinRoot`+`mapBaseChange_injective_adjoinRoot_
      direct`+`pushoutKaehlerSplitStepOption`を2回連鎖させ、
      `V_0=R → V_1=AdjoinRoot g1 → V_2=AdjoinRoot g2`という
      **`d+1`個の単一生成拡大の同時添加**(Faltings の「典型例」
      段落と一致する構成)が実際に組み立てられることを確認した。
      副産物として、2段目の添字族`F`は`e1`(1段目の結果)の型から
      `_`で自動推論させれば明示的に書き下す必要が無いことも分かった
      ——これは`n`段への一般化(`F`をn段分手で書き下す代わりに、
      各段で`_`推論に任せる)の実装コストを大きく下げる発見でもある。
      これで3aで構想した「`V_{n+1}`を`d+1`個の単一生成拡大の反復
      pushoutとして構成する」という道が、**具体例のレベルで完全に
      実証された**——残るのは`n`段(または`d+1`段)の一般化のみ。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **決定的な発見**: Brinon-Conrad Exercise 13.7.4(本ファイル
      2026-09-04の記録で存在は把握していたが、精読していなかった)
      の全文(6ステップ)を実際に読んだところ、Theorem 1.2 の証明の
      核心部分について**上記の`pushoutKaehlerSplitStepOption`(明示的
      なAdjoinRootの同時添加構成)よりもっと直接的な経路**が判明した:

      - **step (1)**: `Ω¹_{B1/A}`が`d+1`個の元で生成される、という
        事実は「second fundamental exact sequence + Nakayama」
        (原文のヒント)で出す——**`B1`を明示的に`d+1`個の生成元の
        AdjoinRoot反復テンソル積として構成する必要は無い**。
        `dim_{k_B}(Ω¹_{k_B/k_A}) = d`(剰余体の`p`-基底の次元、仮定)
        から Nakayama で直接出るので、**どんな有限拡大`B/A`にも
        機械的に適用できる一般命題**になる(構成に依存しない分、
        むしろ`pushoutKaehlerSplitStepOption`より本質的)。
      - **step (2)**: `length_B(Ω¹_{B1/A}) = length_B(B/p^{δ_{B1/A}}B)`
        ——これは Lemma 1.1 そのもの(`falt1CokernelLengthEq`、
        既に**完成・Found**)と同一。
      - **step (4)**: `0→Bₙ₊₁⊗Ω¹_{Bₙ/Aₙ}→Ω¹_{Bₙ₊₁/Aₙ}→Ω¹_{Bₙ₊₁/Aₙ₊₁}→0`
        (完全)+ (1) + elementary divisor theorem から`ker(b)⊇ker(p
        倍写像)`、`length(ker(b∘a))≥length(Bₙ₊₁/p^{βₙ}Bₙ₊₁)`
        (`βₙ=min(1,δₙ/(d+1))`)。左半分の完全性は`KaehlerDifferential.
        exact_mapBaseChange_map`(既存、`falt1_kaehler_length_exact_
        wn1`等で使用済み)+ Lemma 1.1型の単射性条件(`mapBaseChange_
        injective_adjoinRoot_direct`等、既存)でカバーできる見込み。
      - **step (5)**: 判別式の定義から`p^{δₙ-δₙ₊₁}·(Bₙ⊗ₐₙAₙ₊₁)⊆Bₙ₊₁`
        を示し、`coker(b∘a)`が`p^{δₙ-δₙ₊₁}`で消えることを導く——
        これが2026-09-04に記録した「discriminant の塔の乗法性」
        (`Algebra.discr`系、未調査)に相当する、依然として**独立の
        古典的整数論**を要する箇所。
      - **step (6)**: (4)(5)から`δₙ-δₙ₊₁ ≥ βₙ-(d+1)(δₙ-δₙ₊₁)`——
        これは`delta_tendsto_zero`(既に**完成**)の`hrec`仮定と
        厳密に一致する形。

      **`falt1_kaehler_spanFinrank_le`を追加した**(commit分は次項)
      ——step (1) の Nakayama 部分だけを抽出した抽象補題
      (`IsLocalRing.spanFinrank_eq_finrank_quotient`を`Ω[B/A]`に
      適用しただけ)。`Ω[B/A]⊗k_B`の次元が`d+1`以下という**仮定**から
      `Ω[B/A]`が`d+1`個以下で生成されることを結論する——`dim_{k_B}
      (Ω¹_{k_B/k_A})=d`から`Ω[B/A]⊗k_B`の次元が`d+1`以下であることを
      導く部分(局所環のcotangent space理論、`IsLocalRing.
      spanFinrank_maximalIdeal_eq_finrank_cotangentSpace_of_fg`が
      `+1`の出処になる見込み——DVRの`𝔪`は単項なので`spanFinrank 𝔪=1`)
      は次回への課題として残した。

      ★結論(更新): Theorem 1.2 の核心的な困難は、**依然として
      3つの独立した深い部分**(1: Nakayama連鎖の完成、4: 完全列の
      具体的な接続、5: 判別式の塔の乗法性)に分解できることが
      判明した——うち(4)は既存の道具でほぼカバーできる見込みが
      立ち、(1)は今回Nakayama部分だけ切り出せた。(5)(discriminant
      の塔)が依然として**唯一手つかずの独立した古典的整数論**として
      残っている——次回はここ(`Algebra.discr`系のmathlib資産調査)
      から着手するのが最も効率的と判断する。

      ★実際に`Algebra.discr`系を`.cache/mathlib-index.txt`で調査した
      (2026-09-05):
      - `NumberField.isCoprime_differentIdeal_of_isCoprime_discr`
        (`NumberTheory/NumberField/Discriminant/Different.lean`)——
        まさに欲しかった「判別式が互いに素なら differentIdeal も
        互いに素」という事実そのもの。
      - `NumberField.natAbs_discr_eq_natAbs_discr_pow_mul_natAbs_
        discr_pow`——線形独立な拡大の合成体の判別式が各因子の
        判別式の冪の積になる、という塔の乗法性そのもの。

      ★★ただし**両方とも`NumberField`(ℚ上の大域体)専用**に
      定式化されており、Faltings の Theorem 1.2 が要る**局所**の
      設定(完全離散付値環、混標数`(0,p)`)にはそのままでは使えない
      ——mathlib の discriminant-differentIdeal 連携が局所版まで
      カバーしているかは未確認(次回の調査対象)。局所版が無ければ、
      これら大域版の証明技法(`differentIdeal`の`absNorm`との関係
      `NumberField.absNorm_differentIdeal`経由)を局所の Dedekind
      次元1の設定へ移植する、という独立した形式化作業になる見込み
      ——これが Theorem 1.2 の残る核心的な困難の**具体的な次の一手**
      として確定した。

      ★続けて **step (4)(完全列)を完成させた**(`falt1_kaehler_
      length_tower_exact`、commit分は次項)——`Aₙ→Bₙ→Bₙ₊₁`の塔に
      対する`0→Bₙ₊₁⊗Ω[Bₙ/Aₙ]→Ω[Bₙ₊₁/Aₙ]→Ω[Bₙ₊₁/Bₙ]→0`の完全性から
      直ちに長さの加法公式`length(Ω[Bₙ₊₁/Aₙ]) = length(Bₙ₊₁⊗Ω[Bₙ/Aₙ])
      + length(Ω[Bₙ₊₁/Bₙ])`を得た。鍵は`Module.length_eq_add_of_exact`
      (mathlib既存、探索して発見)——右半分の完全性(`KaehlerDifferential.
      exact_mapBaseChange_map`)と全射性(`KaehlerDifferential.
      map_surjective`)は**無条件**(mathlib既存)、`hinj`(左半分の
      単射性)だけがFalt1固有の仮定として残る形にきれいに閉じた。

      これで Exercise 13.7.4 の6ステップのうち、**(2)(=Lemma 1.1)・
      (4)は完成、(1)は部分完成(Nakayama部分のみ)**——残るのは
      (1)の残り(`dim_{k_B}(Ω¹_{k_B/k_A})=d`から`d+1`を導く部分)・
      (3)(elementary divisor theoremでの`ker(b)⊇ker(p倍)`の導出)・
      (5)(discriminantの塔、局所版が無い)の3つ。当初「独立した深い
      部分が3つ」と評価していたが、実際に手を動かしてみると
      **(4)は数行で閉じ**、`Module.length_eq_add_of_exact`のような
      「探せば既にmathlibにある」道具が複数見つかった——引き続き
      (1)(3)(5)の残りを1つずつ潰していく方針が有効と判断する。

      ★続けて step (1) の「`+1`」の出処を特定した(`falt1_dvr_
      cotangentSpace_finrank_eq_one`、commit分は次項)——`Discrete
      ValuationRing`は`IsDiscreteValuationRing.TFAE`(mathlib既存)の
      同値な特徴付けの1つとして**既に`finrank(CotangentSpace)=1`が
      載っていた**(探すまで気付いていなかった、これも「探せば
      mathlibにある」の一例)。

      **ただし残る接続で新しい障害に遭遇した(未解決)**:
      `KaehlerDifferential.exact_kerCotangentToTensor_mapBaseChange`
      から`(maximalIdeal B).Cotangent → Ω[B/An]⊗k_B → Ω[k_B/An]`の
      完全性は得られるが、`kerCotangentToTensor`が`B`-線形・
      `mapBaseChange`が`k_B`(=`ResidueField B`)-線形であるため、
      両者の像・核を`Module.finrank`で比較しようとすると
      `Submodule B _`対`Submodule(ResidueField B)_`という**scalar
      restrictionを跨いだ比較**が必要になる——単純な`rank-nullity`
      の適用では閉じない。これは Exercise 13.7.4 の議論が暗黙に
      使っている「`k_B`上のベクトル空間として見た次元」と「完全列の
      各項が実際に持つ`B`-加群構造」の食い違いを橋渡しする、もう1段
      の配線が要ることを意味する——次回はここ(`Submodule.
      restrictScalars`越しの`finrank`比較、または`LinearMap.
      finrank_range_add_finrank_ker`を`k_B`-線形版に作り直す)から
      再開すること。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **上記の接続を3セッション分の集中的な試行で追い詰めたが、
      まだ閉じていない**——ここまでで判明した正確な障害点を記録する
      (次回、同じ失敗を繰り返さないために):

      - `LinearMap.extendScalarsOfSurjective`(`B→k_B`で`B`-線形写像を
        `k_B`-線形へ持ち上げる、mathlib既存)自体は正しい道具で、
        `kerCotangentToTensor An B k_B : (RingHom.ker(algebraMap B
        k_B)).Cotangent →ₗ[B] k_B⊗[B]Ω[B/An]`を`k_B`-線形に持ち上げる
        ところまでは単独では成功する。
      - **真の障害**: `RingHom.ker(algebraMap B k_B) = maximalIdeal B`
        (`IsLocalRing.ker_residue`から従う)は`rfl`では**ない**
        (実測確認済み)——`algebraMap B (ResidueField B)`は
        `IsLocalRing.residue B`と(`rfl`で)一致するが、`ResidueField
        B`自体は`B⧸maximalIdeal B`と定義的に等しいにもかかわらず、
        この合成全体が`RingHom.ker`を経由すると`rfl`にならない。
      - この命題的な等式を使って`kerCotangentToTensor`の**型**
        (`(RingHom.ker...).Cotangent`)を`CotangentSpace B`(=
        `(maximalIdeal B).Cotangent`)へ書き換えようとすると:
        - `rw [...] at h0`単独(`h0`だけ)は成功するが、
        - **同じ書き換えを`hex0`(`Function.Exact ⇑h0 ⇑g`の形の
          仮定)に対して行おうとすると`"motive is not type correct"`
          で失敗する**(`Function.Exact`の第1引数の型が依存的に
          `algebraMap`に依存するため、`rw`が安全に一般化できない)。
        - `simp only [...] at h0 hex0`も試したが`"simp made no
          progress"`(`hex0`の型のどこに書き換え対象があるか、simp
          の輻輳マッチングでも見つけられない)。
        - `▸`による`LinearEquiv`の手動構成も試したが、motiveの推論に
          同種の理由で失敗した。

      **今後の方針(次回への具体的な引き継ぎ)**: この手の「依存型を
      跨ぐ書き換え」は`rw`/`simp`で押し通すのではなく、
      `Eq.mpr`/`HEq`を明示的に経由するか、あるいは
      `adjoinRootTensorEquivFwd`のときと同じ戦略(mathlibの既存の
      補題を無理に再形成しようとせず、**このセッション固有の設定
      向けに`kerCotangentToTensor`相当の写像を最初から`maximalIdeal
      B`を使って自前で構成し直す**)を取るのが最も見込みが高いと
      判断する——`KaehlerDifferential.kerCotangentToTensor`の証明
      本体(`RingTheory/Kaehler/Basic.lean:775`付近)を読み、
      `maximalIdeal B`版を直接組み立てる方針で次回は開始すること。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **上記の障害を完全に突破し、Exercise 13.7.4 step (1) の本体を
      完成させた**(`falt1_kaehler_finrank_tensor_residueField_le`、
      commit分は次項)。決め手は`rw`/`simp`で`RingHom.ker(algebraMap
      B k_B)`と`maximalIdeal B`を**型として同一視しようとするのを
      完全にやめた**こと——代わりに`Ideal.mapCotangent`(`AlgHom.id`+
      2つの`≤`包含から作れる、mathlib既存の一般構成)を使って**両
      方向の`R`-線形写像**を明示的に構成し、`Ideal.toCotangent`の
      全射性(mathlib既存)を経由して互いに逆写像であることを示す
      ことで、**書き換えを一切経由せず**橋渡しできた
      (`falt1_mapCotangent_maximalIdeal_bijective`)。

      これで:
      1. `falt1_mapCotangent_maximalIdeal_bijective`: `RingHom.ker
         (algebraMap B k_B)`と`maximalIdeal B`(=`CotangentSpace B`)
         を繋ぐ全単射(`R`-線形)。
      2. `LinearEquiv.ofBijective`でこれを`LinearEquiv`化し、
         `kerCotangentToTensor`と合成 →`CotangentSpace B →ₗ[B] Ω[B/An]
         ⊗k_B`という`B`-線形写像`h0'`を得る。
      3. `LinearMap.extendScalarsOfSurjective`(`IsLocalRing.
         residue_surjective`)で`h0'`を`k_B`-線形の`f'`へ持ち上げる。
      4. `Set.range`の**関数レベルの比較**(`Function.Surjective.
         range_comp`、`Submodule`型の違いを一切経由しない)で`f'`の
         像が元の`kerCotangentToTensor`の像と一致することを示し、
         `KaehlerDifferential.exact_kerCotangentToTensor_
         mapBaseChange`(mathlib既存)から`Function.Exact f' g`
         (`g:=mapBaseChange`)を`k_B`線形のまま直接導出。
      5. rank-nullity(`LinearMap.finrank_range_add_finrank_ker`)+
         `falt1_dvr_cotangentSpace_finrank_eq_one`(核の次元≤1)+
         `Submodule.finrank_le`(像の次元≤`d`)を`omega`で結ぶ。

      **教訓(次回以降、同種の「2つの依存的に絡んだ型を同一視したい」
      場面で直接使える)**: `rw`/`simp`で型そのものを書き換えようと
      せず、**両方向の明示的な写像(全単射)を構成してから
      `LinearEquiv.ofBijective`で正式な同型にする**方が、依存型を
      跨ぐ`Function.Exact`のような文脈では遥かに頑健——`Set.range`
      のような**関数レベルの事実**(`Submodule`の型注釈を経由しない)
      で比較を閉じるのも鍵だった。

      これで Exercise 13.7.4 の6ステップのうち、**(1)(=Nakayama
      部分の本体)・(2)(=Lemma 1.1)・(4)が完成**——残るは(1)の
      「`Ω[B/An]⊗k_B`の次元→`Ω[B/An]`が`d+1`個で生成される」への
      最終接続(`falt1_kaehler_spanFinrank_le`の仮定の形と`TensorProduct.
      quotTensorEquivQuotSMul`越しの橋渡しがまだ残る、小さな仕上げ)・
      (3)(elementary divisor theorem)・(5)(discriminantの塔)の3つ。

      ★この最終接続を同じセッションで試みたところ、`TensorProduct.
      quotTensorEquivQuotSMul`(`M`に対し`TensorProduct R(R⧸I)M ≃ₗ[R]
      M⧸I•⊤`、mathlib既存)自体は見つかったが、(a)これも`R`-線形
      (`k_B`-線形ではない、同じ`extendScalarsOfSurjective`パターンで
      対処できる見込み)、(b)`falt1_kaehler_spanFinrank_le`が
      `N:=(⊤:Submodule B Ω[B/An])`という形で立てられており、
      `↥⊤`と`Ω[B/An]`自体の同一視(`Submodule.topEquiv`)がもう1段
      必要、という2つの小さな配線が残っており、`IsScalarTower B k_B
      (TensorProduct B k_B ↥⊤)`のinstance探索でも小さくつまずいた
      ——`falt1_kaehler_finrank_tensor_residueField_le`本体ほどの
      深さの障害ではなく、次回落ち着いて配線すれば閉じる見込みが高い
      (`falt1_kaehler_spanFinrank_le`を`N:=⊤`ではなく`Ω[B/An]`
      直接の形に作り直す方が近道かもしれない、という選択肢も含めて
      検討すること)。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **この最終接続を同じセッションで完成させ、Exercise 13.7.4 step
      (1) が全体として閉じた**(`falt1_kaehler_spanFinrank_le'`・
      `falt1_kaehler_generatedBy_dplus1`、commit分は次項)。上で予告
      した2つの小さな配線(`quotTensorEquivQuotSMul`の`k_B`線形化・
      `Submodule.topEquiv`越しの同一視)を実際に組み立てたところ、
      さらに2つの細かい instance の壁(`Module.IsTorsionBySet.module`
      が自動導出されないため`letI`が要ること、その`IsScalarTower`
      companionが`Module.IsTorsionBySet.isScalarTower`という別名の
      補題として存在すること)に当たったが、いずれも1-2行の`letI`
      追加で解消した。最後に`TensorProduct.quotTensorEquivQuotSMul`
      の型注釈で明示的に`k_B`(`B⧸maximalIdeal B`の代わりに)を
      指定することで`Module kB (TensorProduct B kB Ω[B/A])`の
      instance探索の失敗も解消した。

      **`where`節での結合は`(deterministic) timeout at isDefEq`
      (maxHeartbeats 800000到達)を起こした**——2つの定理を`where`
      で1つの巨大な項にまとめようとすると、依存的なメタ変数の絡み
      合いでelaborationが爆発する(このセッションで既に何度も見た
      パターン、`Falt1Wn1`の名前付きdefパターンと同種の教訓)。
      **2つの独立したtop-level定理として書き、単純な関数適用で
      結ぶだけ**にすると1秒未満で通った——`where`/巨大な結合定理を
      避け、小さな定理を組み合わせるという、このセッション全体を
      通じて繰り返し確認されたLean高速化の鉄則が、ここでも再確認
      された。

      これで Exercise 13.7.4 の6ステップのうち **(1)(2)(4)の3つが
      完全に完成**(以前は(1)が部分完成だったが、今回で本体まで
      到達した)。`falt1_kaehler_generatedBy_dplus1`は
      「`Ω[B/An]`が離散付値環`B`の有限拡大として`d+1`個の元で
      生成される」という、Faltings 原文の該当箇所そのものの
      Lean形式化として単独で成立している。残るは(3)(elementary
      divisor theoremによる`ker(b)⊇ker(p倍)`の導出)・(5)
      (discriminantの塔、局所版)の2つのみ。

      ★★続けて(3)に着手する前に、その**正確な要件**を原文
      (「...which has `(W_{n+1}/pW_{n+1})^{d+1}` as quotient」)まで
      遡って精密に再確認したところ、重要な事実が判明した:
      (3)は`Ω_{V_{n+1}/V_n}⊗W_{n+1}`が**単に**`d+1`個以下で生成
      される(=`falt1_kaehler_generatedBy_dplus1`が与える事実)
      だけでは閉じない——「`(W_{n+1}/pW_{n+1})^{d+1}`への**全射**を
      持つ」という、**`V_n`塔の具体的構成(`d+1`個の同時添加)に
      固有の追加事実**が必要になる(生成元の個数が`≤d+1`であること
      と、剰余体上の次元がちょうど`d+1`であることの両方が揃って
      初めて、直和分解の各因子`W_{n+1}/p^{α_i}`のすべてで`α_i≥1`が
      出る——単なる抽象的なDVR上の加群論だけでは足りない)。

      これは、(3)が**抽象的なDVR一般論として独立に閉じる部分**
      (`falt1_kaehler_generatedBy_dplus1`まで)と、**具体的な`V_n`
      塔の構成(`pushoutKaehlerSplitStepOption`系)に固有の部分**
      (全射性の確認)とに、当初の想定より明確に分かれることを
      意味する——後者は、このセッションで既に構築した
      `pushoutKaehlerSplitStepOption`(`d+1`個の同時添加そのものを
      構成する)と組み合わせて初めて出せる見込みが高い。次回は
      `pushoutKaehlerSplitStepOption`の出力(`Ω[V_{n+1}/V_n}`の
      `d+1`因子分解の**存在**)から、この全射性(`Ω_{Vn+1/Vn}⊗Wn+1`
      が`(Wn+1/p)^{d+1}`へ全射する)を実際に導出することから
      再開するのが筋が良い——`falt1_kaehler_generatedBy_dplus1`
      (抽象DVR一般論)と`pushoutKaehlerSplitStepOption`(具体的な
      塔の構成)という、このセッションで構築した2つの独立した
      道具を**ここで初めて統合する**局面になる見込み。

      ★続けて、この全射性の**単一因子版**を実際に構築した
      (`falt1_omegaAdjoinRoot_surjective_quotient_p`、commit分は
      次項)。Faltings の「典型例」に現れる各因子は`X^p - C a`
      という形の Eisenstein 型多項式——この形の多項式は
      `derivative(X^p-Ca) = C(p)*X^(p-1)`が**常に`p`の倍数**になる
      という一般的な性質を持つ。これと`omegaAdjoinRootQuot`
      (`Ω[AdjoinRoot g/R] ≅ AdjoinRoot g/(derivative)`、既存)を
      組み合わせるだけで、`Ω[AdjoinRoot(X^p-Ca)/R]`が
      `AdjoinRoot(X^p-Ca)/(p)`へ**全射する**ことが直接示せた
      (`Ideal.span_singleton_le_span_singleton`による`(p·x^{p-1})⊆
      (p)`の確認+`Ideal.Quotient.factor`/`factor_surjective`、
      mathlib既存)。

      これで「`d+1`個の同時添加の**各因子**が`(W/p)`へ全射する」
      という、step (3) の全射性要件の**単一因子版**は完成した——
      残るのは、これを`pushoutKaehlerSplitStepOption`(または
      `d+1`因子分解の直和構造)と組み合わせ、「**すべての**因子が
      同時に`(W/p)^{d+1}`へ全射する」という`d+1`因子版へ拡張する
      ことのみ(直和分解された各成分が独立に全射することから、
      積全体の全射性は形式的に従う見込み——`Pi.map`/`Function.
      Surjective.piMap`のような一般論で閉じる可能性が高い)。

      ★確認: `Function.Surjective.piMap`(`Logic/Function/Basic.lean:
      600`)は**実際にmathlibに既存**(`∀i,Surjective(f i) → Surjective
      (Pi.map f)`)——自作する必要は無く、そのまま使える。これで
      「`d+1`因子版への拡張」に必要な一般論はすべて揃っている
      ことを確認した(残るのは`pushoutKaehlerSplitStepOption`の
      `Option ι`添字族と`falt1_omegaAdjoinRoot_surjective_
      quotient_p`(単一因子)を実際に組み合わせる**配線作業**のみ
      ——具体的な`V_n`塔の生成元多項式(`X_i^p-T_i`型`d`個+
      `1`の`p^{n+1}`乗根型1個)を選んで初めて実行できる、次回への
      課題として残す)。

      ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05、
      **step (5)(discriminantの塔、局所版が無い)への潜在的な迂回路
      を発見した**(未着手・戦略的な仮説の段階、次回検討すること)。

      `conductor_mul_differentIdeal`(mathlib既存、`RingTheory/
      DedekindDomain/Different.lean:634`)の正確なステートメント
      ——`(x:B)(hx:Algebra.adjoin K{algebraMap B L x}=⊤):
      conductor(A,x)*differentIdeal(A,B) = span{aeval x(derivative
      (minpoly A x))}`——を再確認したところ、これは`x`が**A[x]がB
      (=積分閉包)へ収束する任意の単生成元**でありさえすれば成り立つ
      **完全に一般的な事実**であり、`falt1_cancelConductorDelta_
      assembled`で使った際に「`x`が`V1`自身の生成元と一致する、
      という本セッションの構成固有の単純化かもしれない」と記録して
      いた懸念は**杞憂だった**——`conductor_mul_differentIdeal`
      自体はどんな単生成元`x`にも適用できる。

      これは、Faltings の「`d+1`個の**同時**添加」ではなく、
      `pushoutKaehlerSplitStepOption`による「**1個ずつ逐次**添加」
      (`d+1`回の微小段)で`V_{n+1}`を構成すれば、**各微小段**で
      `conductor_mul_differentIdeal`(既存・実測済み)がそのまま
      適用できる、という可能性を示唆する——もしこれが機能すれば、
      discriminant の塔という**新しい独立の理論を一切構築せずに**、
      Lemma 1.1 で既に確立した conductor/differentIdeal の道具**だけ**
      でTheorem 1.2の核心を閉じられる見込みがある。

      ★★ただし、これは`delta_tendsto_zero`(既に完成)が要求する
      漸化式`δ_{n+1}≤δ_n-min(1,δ_n/(d+1))/(d+2)`(`d+1`個**同時**の
      構成に固有の形)をそのまま使えないことを意味する——「1個ずつ
      逐次」の構成では、`n`個のマクロ段の代わりに`(d+1)n`個の
      ミクロ段の漸化式になり、**新しい収束の議論を再導出する**
      必要がある(既存の`delta_tendsto_zero`を流用できない可能性が
      高い)。これは discriminant の塔を避けられる代わりに**別の
      新しい数学的議論**(逐次構成での収束)を要する、という
      トレードオフであり、どちらが実際に軽いかは未検証——次回は
      まずこの「逐次構成での収束」が本当に成り立つか(定性的に:
      `d+1`回の逐次添加でも同じ`δ→0`の結論に到達できるか)を、
      小さな具体例か概算で検討することから始めるのが良い。

      ★続けて`cancel_conductor_delta`(既存)の**完全に一般な形**
      (`conductor(Wn,x)*differentIdeal(V1,Wn1) = Ideal.map(...)
      (differentIdeal V0 Wn)`、`falt1_cancelConductorDelta_assembled`
      で使った「`conductor=⊤`」という退化ケースは**この一般形の特殊
      例に過ぎない**、と確認)を精読し、上記の「逐次構成」戦略を
      1段深く分析した。結論: **discriminant の塔を完全に回避できる
      わけではない**——`conductor(Wn,x)`の**サイズ**(=`length(Wn1/
      conductor(Wn,x))`)を`δn`と関係づける下界評価が依然として
      必要で、これが Faltings 原文の「`β`の評価」に相当する部分
      (回避できなかった核心)。

      ただし**朗報**もある: `conductor_mul_differentIdeal`の右辺
      `span{aeval x(deriv(minpoly))}`は、「典型例」の生成元
      (`X^p-T`型)に対しては**具体的に計算可能**(`deriv = p·X^{p-1}`
      、Lemma 1.1 の`differentIdeal_eq_span_derivative`と全く同じ
      パターン、このセッションで何度も使用済み)。つまり必要な
      下界評価は、**外部の(大域体専用の)discriminant理論を移植する
      のではなく**、この具体的な微分の計算+`differentIdeal(V0,V1)`
      自身の計算(同じくEisenstein型多項式の微分から、既存の道具
      だけで閉じる)を組み合わせるだけで**閉じる可能性が高い**——
      これは「新しい独立の古典的整数論」ではなく、「このセッション
      で既に確立した道具(`conductor_mul_differentIdeal`・
      `differentIdeal_eq_span_derivative`・`falt1_length_quotient_
      mul_of_ne_zero`)を、退化ケース(`conductor=⊤`)ではなく
      **一般ケースで具体的に計算し切る**という、範囲は大きいが
      質的には既知の道具の範囲内の作業**に帰着できる見込みが立った。

      次回はこの「`X^p-T`型生成元に対する`conductor(Wn,x)`の具体的な
      下界評価」を、小さな具体例(例えば`d=0`、`p`固定)で実際に
      計算してみることから始めるのが最も見込みが高い——これが
      成功すれば、Theorem 1.2 の残る核心的な困難(従来「discriminant
      の塔」と呼んでいたもの)が、**このセッションの既存の道具の
      範囲内で解決できる**ことが確定する。

      ★★step (3) の「`d+1`因子版への拡張」に要る**残り**の一般論
      (単一因子の全射性を`B`へ base change する部分)も確認した:
      `LinearMap.rTensor_surjective`/`LinearMap.lTensor_surjective`
      (mathlib既存、`f`が全射なら`f⊗id`も全射)——`Function.
      Surjective.piMap`と合わせ、step (3) の`d+1`因子版に**必要な
      一般論はすべてmathlibに既存**であることが確定した。残るのは
      `pushoutKaehlerSplitStepOption`の`Option ι`出力・`falt1_
      omegaAdjoinRoot_surjective_quotient_p`(単一因子)・これら3つの
      一般論(`piMap`・`rTensor_surjective`・`lTensor_surjective`)を
      実際に組み合わせる**配線作業のみ**——具体的な`V_n`塔(`X_i^p-T_i`
      型`d`個+冪根型1個の生成元)を選んで初めて実行できる。

      ★★★今セッションの到達点のまとめ(次回への引き継ぎ):
      Theorem 1.2 の証明は Brinon-Conrad Exercise 13.7.4 の6ステップ
      に分解でき、うち**(1)(2)(4)は完全に完成**、**(3)は単一因子版が
      完成し、d+1因子版に要る一般論もすべて確認済み**(残るのは
      具体的な`V_n`塔を選んでの配線のみ)、**(5)は「discriminantの塔」
      という独立理論の移植ではなく、既存の`conductor_mul_
      differentIdeal`・`differentIdeal_eq_span_derivative`による
      **具体的な計算**に帰着する見込みが立った(まだ実行していない)。
      残る2つの困難は、いずれも「新しい数学理論の輸入」ではなく
      「具体的な`V_n`塔を選んでの計算・配線」という**質的に同じ
      種類の作業**に帰着したことが、今セッション最大の戦略的成果。

      ★★★★上記の「discriminant を回避できるかもしれない」という
      期待(2026-09-05の「逐次構成」の記録)を、`cancel_conductor_
      delta`の一般形をさらに1段掘り下げて検証したところ、**過度に
      楽観的だった**ことが判明した——正直に訂正する。`condWnx *
      Jn = spanDeriv`(`conductor_mul_differentIdeal`)で`Jn :=
      differentIdeal(Wn,Wn1)`を計算するには、`x`の**`Wn`上の**
      最小多項式が要る。`falt1_cancelConductorDelta_assembled`が
      「`conductor=⊤`」という退化結果になったのは、`x`が`V1`の
      生成元と**一致する**ため、`Wn`上の最小多項式も**同じ`X^n-π`**
      になり、`Wn[x]=Wn1`が(ズレ無く)ちょうど成り立ってしまう
      からだった。**真に異なる生成元**を使えば退化は避けられるが、
      その場合`x`の`Wn`上の最小多項式(`=`「`X^n-π`が`Wn`上でどう
      分解・分岐するか」)は**それ自体が独立の計算対象**になり、
      これはまさに古典的な塔の分岐理論(discriminant/conductor理論)
      **そのもの**——結局、名前を変えただけで**同じ核心的困難に
      戻ってくる**ことが分かった。「既存道具の範囲内の計算」という
      評価も訂正が要る: 既存道具(`conductor_mul_differentIdeal`等)
      は計算の**枠組み**を提供するが、**入力となる`x`の`Wn`上の
      最小多項式(またはそれに相当するramification情報)は依然として
      独立に確立する必要がある**——これが Theorem 1.2 の核心の
      困難の**正確な所在**である、という当初(2026-09-04)の評価に
      結局回帰した。次回はこの認識(discriminant/ramification理論の
      核心は回避できない)を前提に、それでも**局所の設定でどこまで
      計算可能か**(`Algebra.discr`の局所版が無いなら、`differentIdeal`
      の塔公式を繰り返し適用しながら`x`の最小多項式の次数・分岐指数
      を直接追跡する、という原始的だが着実な道)を検討することから
      始めるのが誠実だと判断する。

      ★念のため`ramificationIdx`と`differentIdeal`(または`discr`)を
      結ぶ定量的な補題(古典的な「`e-1 ≤ v(𝔡) ≤ e-1+v(e)`」型の
      different指数の公式)が mathlib に無いか`.cache/mathlib-index.
      txt`で最終確認したが、`ramificationIdx`に言及する92件の中に
      `differentIdeal`/`discr`との連携は1件も見当たらなかった——
      これでmathlibにこの定量的な橋渡しが無いことを複数の角度から
      確認済み(discriminant系・conductor系・ramificationIdx系の
      いずれからも同じ結論)。Theorem 1.2 の残る核心は、確かに
      「mathlibに無い独立した古典的整数論」であることが最終確認
      できた——今回の探索は空振りではなく、**探すべき場所を
      使い切った**という意味で意義がある。

      ★step (3) の`d+1`因子版への配線も同じセッションで試みたところ、
      **単一因子の全射性を`B`へ base change する**部分で新しい技術的な
      詰まりに当たった(未解決、次回への課題として記録):
      `falt1_omegaAdjoinRoot_surjective_quotient_p`は`≃+`
      (`omegaAdjoinRootQuot`経由)ベースなので**加法群としての全射**
      にすぎず、`LinearMap.lTensor`(`AdjoinRoot g`-線形写像を要求)へ
      渡すには`AdjoinRoot g`-線形であることを別途示す必要がある——
      `omegaAdjoinRootQuot`が実際に`AdjoinRoot g`-線形かどうかは
      `exact?`では自動確認できず(未証明、おそらく真だが`omega_
      quotient_eq_derivative_span`の構成まで遡って確認する必要が
      ある)。次回はここ(`omegaAdjoinRootQuot`を`≃+`ではなく
      `≃ₗ[AdjoinRoot g]`として再構成する、または全射性の証明を
      最初から`LinearMap`ベースで書き直す)から再開すること。

      (この段落で構想した代替路は上で実際に`falt1_differentIdeal_
      tower_length`として確立・commit済み——詳細は上記参照。project内
      の`differentIdeal_tower_diamond`は同じmathlib補題を2回使う
      「菱形」版であって、対応物というより`falt1_differentIdeal_
      tower_length`の親戚にあたる、との理解で確認済み。)

      ★★2026-09-05、上記の「次回はここから」を実際に再開し、**解決
      した**。`omegaAdjoinRootQuot`(`≃+`のみ)を直接`AdjoinRoot g`-線形
      と示す代わりに、既に`≃ₗ[AdjoinRoot g]`である`omega_quotient_eq_
      derivative_span`を土台に**同じ定理を線形版として再構築**する
      方針に切り替えた:
      - `falt1_omegaAdjoinRoot_surjective_quotient_p_lin`:
        `Ω[AdjoinRoot g⁄R] →ₗ[AdjoinRoot g] AdjoinRoot g⧸(p)`への全射。
        `Submodule.mapQ`+`Submodule.factor_surjective`で`Ideal.
        Quotient.factor`の線形版を構成(`hle`の証明は`≃+`版と同一の
        微分計算)。
      - `falt1_omegaAdjoinRoot_surjective_quotient_p_baseChange`:
        上記を`LinearMap.baseChange B`で`B`へ base change し、
        `TensorProduct (AdjoinRoot g) B Ω[AdjoinRoot g⁄R] →ₗ[B]
        TensorProduct (AdjoinRoot g) B (AdjoinRoot g⧸(p))`への全射
        (`LinearMap.baseChange_surjective`)を得た。
      - `falt1_omegaAdjoinRoot_surjective_quotient_p_baseChange'`:
        右辺のテンソル積を`TensorProduct.quotTensorEquivQuotSMul`+
        `TensorProduct.comm`(因子の順序が`quotTensorEquivQuotSMul`の
        想定と逆だったため)で`B⧸(I•⊤ : Submodule (AdjoinRoot g) B)`
        (単なる`B`の商加群)に繋ぎ直した。`Module`インスタンスは
        `Module.isTorsionBySet_quotient_ideal_smul`から`letI`でその場
        構成(グローバル登録はしない)。一度`Nonempty (... →+ ...)`
        という誤った目標形で`⟨map, ?_⟩`のanonymous constructorが
        「加法性の証明」枠に誤マッチする(先の`Prod`誤マッチと同種の
        バグ)エラーに当たったが、`∃ ρ, Function.Surjective ρ`という
        既に機能していた形に戻すことで解決した。3定理とも`lean_check`
        で確認後、`KaehlerAux.lean`に追記・`lake build`成功(sorry
        無し)。これで Exercise 13.7.4 step (3) の**単一因子の線形
        base change**が完成し、残るは`pushoutKaehlerSplitStepOption`
        との統合による`d+1`因子版への配線のみ。

      ★★2026-09-05、続けてその「配線」の**帰納の1段**を実際に確立
      した。`pushoutKaehlerSplitStepOption`は`Ω[B/R] ≃ₗ[B] (∀i:Option
      ι,...)`という`Option ι`添字の積への同型を与える——これと(a)
      新規因子`C`についての全射`φC`(`falt1_omegaAdjoinRoot_surjective_
      quotient_p_baseChange`を渡す想定)・(b)前段までの各因子`F i`に
      ついて`B`へ運んだ後の全射`φ i`の**成分ごとの全射性**を
      `LinearMap.piMap`(mathlibの、成分ごとの線形写像から積型の間の
      線形写像を作る道具)+`Function.Surjective.piMap`(純粋な集合論の
      「成分ごと全射→積全体で全射」)で束ね、`Ω[B/R]`から`(∀i:Option
      ι, Option.elim i QC Q)`への全射を得る`falt1_pushoutKaehlerSplit
      StepOption_surjective`を完成させた。副産物として:
      - `falt1_piMap_surjective`: `LinearMap.piMap`の全射性
        (`LinearMap.coe_piMap`で強制関数の記述を`Pi.map`に落とし、
        `Function.Surjective.piMap`を適用するだけ)。
      - `falt1_optionElim_addCommGroup`/`falt1_optionElim_module`:
        `Option.elim i QC Q`という**依存する**族に対する`AddCommGroup`/
        `Module`インスタンス。`match`の各枝で`Option.elim`が定義通り
        簡約されることを型クラス探索は自動では認識しない
        (`inferInstance`だけでは`AddCommGroup (none.elim QC Q)`が
        解決できない)ため、`show`で先に簡約後の型へ変換してから
        `infer_instance`する、という一手間が必要だった——新しい失敗形
        として`tools/lean-idioms.md`に足す価値がある。
      3定理とも`lean_check`で確認後、`KaehlerAux.lean`に追記・`lake
      build`成功(sorry無し)。これで Exercise 13.7.4 step (3) の
      「帰納の1段」が一般形で完成した。残る作業は、この帰納を
      **実際の`V_n`塔**(`d+1`個のEisenstein型生成元を`pushoutKaehler
      SplitStepOption`で`n`回連鎖)に対して`n`回回し、`falt1_omegaAdjoin
      Root_surjective_quotient_p_baseChange`を各段の`φC`として実際に
      渡す配線——`falt1_pushoutKaehlerSplitStepOption_adjoinRoot_
      example`(具体例、`True`で終わる骨格のみ)を土台に、そこへ実際に
      全射性の主張を載せる作業として次回に持ち越す。なお`pushoutKaehler
      SplitBase`が使う`ι=Fin 1`の出発点(`R`自身をダミーの1因子とする)
      は`Ω[R/R]=0`を寄与するだけなので、対応する`QC`/`Q`はダミーの
      自明加群(例えば`PUnit`)を渡せば良いはずで、実質的な障害には
      ならない見込み(未検証)。

      ★2026-09-05、続けてこの「未検証」を実際にREPLで試したところ、
      **`PUnit`自体は問題無く動いた**(`Subsingleton`な`Ω[R/R]`から
      `PUnit`への零写像が全射であることは`⟨0,Subsingleton.elim _ _⟩`
      で即座に閉じる、`PUnit.{u+1}:Type u`という宇宙も明示すれば
      問題無い)が、**別の箇所で本物のinstance diamondに当たった**:
      `falt1_pushoutKaehlerSplitStepOption_surjective`の`φ`が要求する
      型は`(F i).lift (B:=B2)).carrier`上の`Algebra`(`RAlgOver.lift`が
      `((algebraMap B1 B).comp (algebraMap x.carrier B1)).toAlgebra`
      として構成する、合成RingHom経由の instance)だが、`(F i).carrier
      = R`・`B2 = AdjoinRoot(g1.map(algebraMap R R))`という具体例では
      `B2`に**既に`AdjoinRoot.instAlgebra`という自然な`Algebra R B2`
      instanceが存在**しており、両者は`algebraMap R R = RingHom.id R`
      (`Algebra.algebraMap_self`、命題として真だが`rfl`ではない)経由
      でしか繋がらない——`letI`で`.lift.algT`を明示的に登録しても
      `synthesized type class instance is not definitionally equal`
      で弾かれた(`tools/lean-idioms.md` #1「instances透明度で型が
      合わない」の**まさに典型例**が、今回はダミー因子`R`という
      最も単純なケースでさえ発生することを確認)。これは`falt1AdjoinRoot
      Algebra`(`def`として分離、「衝突を避けるため意図的に」と明記
      されている)が対処したのと同種の壁——**本物の直し方**は、
      `algebraMap R R = RingHom.id R`の証明を経由して`φ`を`cast`する
      か、あるいは`Fdummy`を`⟨R⟩`ではなく最初から`AdjoinRoot.instAlgebra`
      と同じ経路で構成される項に取り替える(ダミー因子であっても
      「本物の`AdjoinRoot`」として扱う)ことだと考えられるが、
      いずれも次回への持ち越しとする——ここで無理に押し切って壊れた
      ものを残すより、`falt1_pushoutKaehlerSplitStepOption_surjective`
      という**一般形の帰納の1段**(既に`lake build`で確認済み・commit
      済み)を確実な到達点として確定させる方を優先した。

      ★★★2026-09-05、続けてこのinstance diamondを**実際に解消した**。
      原因は診断通りではなく、もう1段単純だった: `φ`/`hφ`を`have`で
      **独立に型注釈してから**渡すと、その型注釈の elaboration が
      (`letI`で明示登録したにも関わらず)`AdjoinRoot.instAlgebra`を
      独自に見つけてしまい、理論側が期待する`.lift.algT`と食い違う。
      **直し方**: `φ`/`hφ`を独立に`have`/`set`せず、`falt1_pushoutKaehler
      SplitStepOption_surjective`の呼び出しの**引数の位置に無名関数
      として直接書く**(`fun i => 0`・`fun i y => ⟨0, Subsingleton.elim
      _ _⟩`)——こうすると型が**呼び出し先(理論の`φ`引数の型)から
      直接推論**され、独立した型注釈によるinstance探索の分岐が
      そもそも発生しない。`falt1_pushoutKaehlerSplitStepOption_
      adjoinRoot_surjective_example`として`KaehlerAux.lean`に追記・
      `lean_check`確認後`lake build`成功。これで Exercise 13.7.4
      step (3) の帰納の1段が、**実際のAdjoinRoot具体例(1生成元+
      ダミー因子)に対して動くことを end-to-end で確認**できた——
      `d+1`個の実生成元への拡張(`Fdummy`側を本物の前段の族に
      置き換えて`n`回連鎖する)への道筋がこれで裏付けられた。
      新しい教訓として`tools/lean-idioms.md`にも追記した(「独立
      `have`型注釈 vs 呼び出し引数位置での無名関数」という、既存の
      `#1`/`#33`とは別角度の直し方)。

      ★★★2026-09-05、続けて`falt1_kaehler_length_exact_wn1_kernel`の
      docstringが「次はこれに進む」と記録していた**合流**を実際に
      実行した。`falt1_kaehler_length_exact_wn1_cokernel`(右辺第2項=
      `Wₙ₊₁⧸differentIdeal Wₙ Wₙ₊₁`)+`falt1_kaehler_length_exact_wn1_
      kernel`(右辺第1項=`Wₙ₊₁⧸(differentIdeal V0 Wₙ).map(...)`)を
      合わせると、`length Ω¹_{Wₙ₊₁/V0}`が`falt1_differentIdeal_tower_
      length`(differentIdealの塔公式)の右辺と**完全に同じ2項の和**
      になっていることに気づき、`falt1_kaehler_length_exact_wn1_full`
      として3つを貼り合わせた(`add_comm`+`rw`のみ、`lean_check`一発
      でOK)。結論:
      ```
      length Wₙ₊₁ Ω¹_{Wₙ₊₁/V0} = length Wₙ₊₁ (Wₙ₊₁⧸differentIdeal V0 Wₙ₊₁)
      ```
      これは`falt1CokernelLengthEq`(Lemma 1.1)を`V0→Wₙ₊₁`に直接適用
      した場合と**同じ結論**だが、**異なる経路**——`falt1_theorem12_
      kaehler_length`が要求していた「`Wₙ₊₁`の`V0`上の単項生成元`y`
      (一般に存在するか未確認、と記録していた仮定)」が**そもそも
      不要**で、代わりに`Wₙ`の`V0`上の生成元`w`と`Wₙ₊₁`の`Wₙ`上の
      生成元`x`という、**それぞれの帰納的構成が自動的に供給する
      弱い前提**だけで足りる。これで falt1_theorem12_kaehler_length_
      eq_differentIdeal が抱えていた「一般存在性が未確認」という
      条件付き完成の障害を、この経路では回避できることが分かった
      ——`lake build`で確認後、次回はこれを実際の`V_n`/`W_n`帰納
      (`w`・`x`を各段の具体的な生成元に、`hsep`・`hne`を各段で確認)
      に接続し、δ_n の漸化式(`hrec`)の**左側**(`length Ω¹`から
      `differentIdeal`への変換)を完全に埋める作業に進む。残る
      核心の困難(step 5、`differentIdeal Wₙ Wₙ₊₁`の**下からの評価**、
      すなわち`β`項の由来)は、この合流だけでは解消しない——引き続き
      新しい局所分岐理論を要する。

      ★★★★2026-09-05、続けて**step (5)の核心(`β`項の由来)に
      初めて実際に着手した**。`conductor_mul_differentIdeal`+
      `conductor_eq_top_of_adjoin_eq_top`(既存、`differentIdeal_eq_
      span_derivative`として既にwrapされていた)から`differentIdeal
      Wₙ Wₙ₊₁ = span{n·x^(n-1)}`(`x`=Eisenstein型`g=X^n-Cπ'`の根、
      微分の計算は`falt1_hspan_eq`の証明に既に現れていたが、独立の
      定理としては未抽出だった)——ここまでは既存の断片の組み合わせ。
      **新規に構成したのは**、「`x`が`Wₙ₊₁`の一様化元である」
      (Eisenstein拡大が全分岐であることの核心、mathlibに無いと
      確認済みの事実)の**代数的な言い換え**を、mathlibに直接の道具が
      無いため**自前で証明した**こと:
      - `falt1_adjoinRoot_quotient_root_eq_residue`:
        `AdjoinRoot(X^n-Cπ')⧸(root) ≃ₐ[Wₙ] Wₙ⧸(π')`(第一同型定理を
        自前の`AlgHom`(`AdjoinRoot.liftAlgHom`、根を`0`へ送る)に
        適用。核の計算は`p∈ker`なら`p=C(p.coeff 0)+X·q`という初等
        分解と`x^n=π'`(Eisenstein関係式そのもの)を組み合わせて
        `p∈(root)`を導く、という完全に初等的な議論のみで閉じた)。
      - `falt1_adjoinRoot_quotient_root_isField`: 上記の系として、
        `Wₙ`が DVR で`π'`が一様化元なら`Wₙ⧸(π')`が体(剰余体)——
        ゆえに`AdjoinRoot g⧸(root g)`も体、すなわち`(root g)`が
        `AdjoinRoot g`の**極大イデアル**であることを得た
        (`MulEquiv.isField`で体の性質を同型で輸送するだけ)。
      これは「根が生成するイデアルは極大」という、全分岐の**核心
      部分**を初めて実際に証明したもの——`tools/lean-idioms.md`の
      教訓通り`AdjoinRoot.liftAlgHom`の`0`引数で数回宇宙/型推論の
      罠に当たったが(`z0`を`set`で独立に固定することで解決)、
      それ以外は完全に見通し通りに閉じた。`lean_check`確認後、
      `KaehlerAux.lean`に追記。次回: (a) `Wₙ₊₁ ≃ₐ[Wₙ] AdjoinRoot g`
      経由でこの結果を実際の`Wₙ₊₁`へ輸送し、(b) `Wₙ₊₁`が局所環である
      ことと合わせて`(root g)=maximalIdeal Wₙ₊₁`(極大ideal→
      局所環の**唯一の**極大idealなので一致)を導き、(c)
      `IsDiscreteValuationRing.length_quotient_pow_maximalIdeal`で
      `length(Wₙ₊₁⧸(x)^(n-1))=n-1`を得て、`length(Wₙ₊₁⧸differentIdeal
      Wₙ Wₙ₊₁) = length(Wₙ₊₁⧸(n)) + (n-1) ≥ n-1`という、`β`項の
      **具体的な下限**を確立することに進む。

      ★★★★★2026-09-05、続けて(a)(b)を**同じセッションで完成させた**。
      `falt1AdjoinRootEquivIntegralClosure`(既存、`AdjoinRoot g ≃ₐ[Wₙ]
      Falt1Wn1 Wn Wn π' n`)で`falt1_adjoinRoot_quotient_root_isField`を
      実際の`Wₙ₊₁`(`Falt1Wn1 Wn Wn π' n`、既存の名前付き`def`を
      `V0:=Wn`で「自己適用」して再利用——巨大な入れ子式の再記述を
      避けつつ独立の名前付き`def`を新設する必要が無かった)へ輸送し、
      さらに**「唯一性」まで完全に証明した**:
      - `falt1_falt1Wn1_uniqueMaximalIdeal`: 任意の極大イデアル`𝔪`に
        ついて、`Ideal.isMaximal_comap_of_isIntegral_of_isMaximal`
        (既存、整拡大では極大idealの引き戻しは極大)+`Wₙ`が局所環
        であることから`𝔪`の`Wₙ`への引き戻しが`maximalIdeal Wₙ=(π')`
        と一致し、`x^n=π'∈𝔪`(`x`=根の像)・`𝔪`が素なので`x∈𝔪`、
        `(x)`が極大(前段の結果)なので`(x)⊆𝔪⊊Wₙ₊₁`から`(x)=𝔪`、と
        いう初等的な議論で、`(x)`が**唯一の**極大イデアルであることを
        示した。
      - `falt1_falt1Wn1_isLocalRing`: 上記の系として`IsLocalRing.
        of_unique_max_ideal`を適用し、`Falt1Wn1 Wn Wn π' n`(Falt1の
        実際の`Wₙ₊₁`構成)が**局所環**であることを得た。
      これで**「Eisenstein拡大は全分岐である」という、mathlibに
      存在しないと確認済みだった古典的整数論の事実を、完全に自前で
      証明し切った**——`tools/lean-idioms.md`の教訓(`Falt1Wn1 Wn Wn
      π' n`の「自己適用」による巨大な入れ子式回避、`haveI ‹...›`
      による名前付き`def`↔生の`integralClosure`式の instance 橋渡し)
      以外は完全に見通し通りに閉じた。`lean_check`確認後
      `KaehlerAux.lean`に追記。次回: `IsDedekindDomain + IsLocalRing
      ⟹ IsDiscreteValuationRing`(TFAE経由、既に`IsDiscreteValuationRing.
      TFAE`で確認済みの同値性)を適用して`Falt1Wn1 Wn Wn π' n`が実際に
      DVRであることを導き、`IsDiscreteValuationRing.length_quotient_
      pow_maximalIdeal`で`length(Wₙ₊₁⧸(x)^(n-1))=n-1`という`β`項の
      **具体的な下限**を確立する——これがstep (5)の最後の技術的な
      接続になる見込み。

      ★★★★★★2026-09-05、続けてこの「最後の技術的な接続」を
      **同じセッションで完成させた**。当初想定した`IsDiscreteValuation
      Ring.TFAE`経由ではなく、**もっと直接的な経路**があると判明した:
      mathlibの`class IsDiscreteValuationRing`は`extends IsPrincipal
      IdealRing, IsLocalRing`に加え、**もう1つの隠れたフィールド**
      `maximalIdeal ≠ ⊥`(体でないこと)を要求するだけ——`IsPrincipal
      IdealRing`(既存の`falt1_isPrincipalIdealRing_of_finite_ext_of_
      DVR`)+`IsLocalRing`(前段の結果)を`constructor`で組み合わせ、
      残る`maximalIdeal≠⊥`を`x≠0`(`x^n=algebraMap π'≠0`、後者は
      `Module.IsTorsionFree`から「`π'`が正則なら`Falt1Wn1...`上でも
      正則に作用する」`Module.IsTorsionFree.isSMulRegular`経由)から
      示すだけで**TFAEを経由せず直接閉じた**。
      - `falt1_falt1Wn1_isDiscreteValuationRing`: `Falt1Wn1 Wn Wn
        π' n`が`IsDiscreteValuationRing`であることの完全な証明
        (前段の`uniqueMaximalIdeal`系2定理を呼ぶのではなく、`x`の
        一貫性を保つため同じ構成を1つの証明にまとめて自己完結させた)。
      - `falt1_falt1Wn1_length_quotient_maximalIdeal_pow`: 上記の
        DVR instanceを使い、`IsDiscreteValuationRing.length_quotient_
        pow_maximalIdeal`(mathlib既存)を直接適用して`length(Wₙ₊₁⧸
        maximalIdeal^k)=k`(**任意の`k`で**)を得た——`k:=n-1`とすれば
        まさに探していた`β`項の下限の土台。
      これで**「Eisenstein拡大が全分岐である」という事実から、
      具体的な長さの下限(`length(Wₙ₊₁⧸(x)^(n-1))=n-1`)までを、
      完全にmathlib既存の道具(TFAEすら経由せず)で導き切った**——
      Theorem 1.2 step (5)の技術的な核心はこれで確立できたと言える。
      `lean_check`確認後`KaehlerAux.lean`に追記。残る接続作業
      (`differentIdeal Wₙ Wₙ₊₁=span{n·x^(n-1)}`とこの`length(Wₙ₊₁⧸
      maximalIdeal^(n-1))=n-1`を`falt1_length_quotient_mul_of_ne_zero`
      で貼り合わせ`length(Wₙ₊₁⧸differentIdeal Wₙ Wₙ₊₁)=length(Wₙ₊₁⧸(n))
      +(n-1)`という完全な等式にすること、さらにこれを`hrec`の`β`項へ
      正式に接続すること)は次回に持ち越す。

      ★2026-09-05、この「残る接続作業」に着手し、`differentIdeal_eq_
      span_derivative`(既存)+`falt1BaseChangeGeneratorFull`(既存、
      `hxadjoin`・`hxminpoly`込みの生成元`x`を直接供給)+`falt1_
      fieldLevel_adjoin_top_of_ringLevel_minpoly`(既存、体レベルの
      adjoin条件への橋渡し)を組み合わせて`differentIdeal Wₙ Wₙ₊₁=
      span{n·x^(n-1)}`を直接示そうとしたが、`differentIdeal_eq_span_
      derivative`の`hw`引数(体レベルのadjoin条件)が要求する暗黙の
      `[IsIntegralClosure Wₙ₊₁ Wₙ L]`等のinstanceと、`falt1_fieldLevel_
      adjoin_top_of_ringLevel_minpoly`が実際に供給する**具体的な**
      `L:=AdjoinRoot(g.map(FractionRing Wₙ))`に対するinstanceとが
      自動では繋がらず、`isDefEq`が`maxHeartbeats 1000000`でも
      timeoutする(このファイルが`falt1_cancelConductorDelta_
      assembled`等で繰り返し記録してきた「巨大な入れ子式」由来の
      elaboration負荷の**同じ症状**、しかし今回は表示だけでなく
      実際の型検査自体が重い、より深刻なケース)。`falt1_hspan_eq`
      (既存)がまさにこの橋渡しを内部で行っているので、そちらの
      instance登録の手順(`haveI`の並び)を精読して踏襲すれば解決
      できる可能性が高いが、この手順を正確に再現するには
      `falt1_hspan_eq`の証明全体(150行超)をもう一段精読する時間が
      必要——次回はそこから再開する。無理に押し切って壊れたものを
      残すより、`falt1_falt1Wn1_isDiscreteValuationRing`・`falt1_
      falt1Wn1_length_quotient_maximalIdeal_pow`という確実な到達点
      (共にcommit済み)を維持する方を優先し、この接続作業自体は
      コミットしなかった。

      ★★★2026-09-05、続けてこの`isDefEq`timeoutを**実際に解消した**。
      原因は診断通りだった: `differentIdeal_eq_span_derivative`を
      呼ぶ**前**に、`falt1_hspan_eq`が内部で毎回行っていた`Fact
      (Irreducible ...)`・`FiniteDimensional`・`Algebra.IsSeparable`
      の事前登録(`g`側=`Wₙ→Wₙ₊₁`のfield level)を**省略していた**
      ことが原因だった——これらを呼び出し前に明示的に`haveI`すれば、
      暗黙のinstance探索が正しく解決し、`maxHeartbeats 1000000`でも
      timeoutせず**3秒程度**で閉じた(`set_option maxHeartbeats`
      すら実質不要なレベル)。`falt1_differentIdeal_Wn_Wn1_eq_span_
      deriv`として
      ```
      differentIdeal Wₙ Wₙ₊₁ = Ideal.span {n · x^(n-1)}
      ```
      を`lean_check`で確認後`KaehlerAux.lean`に追記した——これで
      step (5)の`β`項の由来(`differentIdeal`の閉じた式)と、
      前々回確立した具体的な長さの下限(`length(Wₙ₊₁⧸maximalIdeal^k)
      =k`)の**両方**が、同じ`x`(`falt1BaseChangeGeneratorFull`が
      供給する生成元)を経由して**言葉が揃った**。残る最後の接続
      (`Ideal.span{n·x^(n-1)}=Ideal.span{n}*Ideal.span{x}^(n-1)`
      への分解+`falt1_length_quotient_mul_of_ne_zero`で`length(Wₙ₊₁⧸
      differentIdeal Wₙ Wₙ₊₁)=length(Wₙ₊₁⧸(n))+(n-1)`という完全な
      等式にまとめること)は次回に持ち越すが、必要な材料は**すべて
      揃った**。

      ★★★★★★★2026-09-05、この「残る最後の接続」を**同じセッション
      で完成させた**。鍵となる発見は、`falt1_falt1Wn1_isDiscreteValuation
      Ring`が独自に構成していた`e2`(`falt1AdjoinRootEquivIntegralClosure`
      経由)を再利用する代わりに、`x`(`falt1BaseChangeGeneratorFull`
      供給)自身の`hxadjoin`・`hxminpoly`から**直接**`AdjoinRoot(minpoly
      Wₙ x) ≃ₐ[Wₙ] Wₙ₊₁`を構成する`adjoinRootMinpolyEquiv`
      (`IsAdjoinRootMonic.mkOfAdjoinEqTop`経由、既存)を使えば、
      「根が`x`に写る」ことが`adjoinRootMinpolyEquiv_root`+
      `IsAdjoinRoot.adjoinRootAlgEquiv_apply_root`(共に既存)から
      **無条件に**従うと気づいたこと——これで`falt1_differentIdeal_
      Wn_Wn1_eq_span_deriv`の`x`(`differentIdeal`の閉じた式の根拠)
      と、全分岐(極大idealの一意性)の議論の`x`を**別々に構成して
      一致を証明する**という難所を、そもそも回避できた。
      `falt1_theorem12_length_differentIdeal_eq_length_quotient_n_add`
      として、`falt1_differentIdeal_Wn_Wn1_eq_span_deriv`(閉じた式)
      +全分岐の議論(自己完結、`falt1_falt1Wn1_isDiscreteValuationRing`
      の議論を`x`の元で再構成)+`Ideal.span_singleton_mul_span_
      singleton`・`Ideal.span_singleton_pow`(mathlib既存、`span{a·b}
      =span{a}*span{b}`・`span{a}^k=span{a^k}`)+`falt1_length_
      quotient_mul_of_ne_zero`+`IsDiscreteValuationRing.length_
      quotient_pow_maximalIdeal`を貼り合わせ、
      ```
      length Wₙ₊₁ (Wₙ₊₁ ⧸ differentIdeal Wₙ Wₙ₊₁)
        = length Wₙ₊₁ (Wₙ₊₁ ⧸ span{n}) + (n - 1)
      ```
      を確立した(`ℕ∞`の引き算の橋渡しに`ENat.coe_sub`が要った、
      小さいが新しい罠)。`lean_check`で約12秒(この規模の証明として
      は高速)で確認後`KaehlerAux.lean`に追記した。これで**Eisenstein
      拡大が全分岐であるという事実から、Theorem 1.2 の漸化式(`δₙ→0`)
      の`β`項の具体的な下限(`n-1`、`n=p`なら`p-1`)までを、完全に
      自前で証明し切った**——step (5)の技術的な核心はこれで**完成**
      した。次回: この結果を実際の`V_n`塔の帰納(`δ_n`・`δ_{n+1}`を
      `differentIdeal V0 Wₙ`・`differentIdeal V0 Wₙ₊₁`の言葉で結ぶ、
      `falt1_kaehler_length_exact_wn1_full`との統合)に接続し、
      `delta_tendsto_zero`が要求する`hrec`の形へ正式に持ち上げる
      作業に進む——これがTheorem 1.2の証明の最後の大きな統合作業。

      ★★★★★★★★2026-09-05、この「次回」に実際に着手し、**重要な、
      過度に楽観的だった見立ての訂正**に至った。原文(本ファイル
      152-173行目に逐語記録済み)を読み直すと:「If we denote the
      different **of W_n over V_n**」——`δ_n`は`Wₙ`と`Vₙ`という
      **2つの独立な塔**の差分であり、`Vₙ`(このセッションで構築した
      Eisenstein型単生成塔がまさにこれに相当する——「典型例」段落の
      `V_{n+1}/V_n`の記述と完全に一致)と、`Wₙ`(**`Vₙ⊗_V W`の
      正規化**、`W`はTheorem 1.2が扱う**任意の**almost étale
      covering)は**別物**。今回確立した`length(Wₙ₊₁⧸differentIdeal
      Wₙ Wₙ₊₁)=...`は、実は**`Vₙ→Vₙ₊₁`塔自身**(Faltings記法での
      `Vₙ`側)についての結果であり、`δₙ=different(Wₙ,Vₙ)`が要求する
      **`Wₙ`側**(`W`の情報を担う、almost étale coveringの正規化)
      は依然として**未着手**——`W=V`(自明なalmost étale covering)
      という**退化した特殊ケース**でなら`δₙ≡0`となり漸化式は
      空虚に真になるが、これは§2-4で確認済みの「非空虚性の欠如」
      そのものであり、G9の精神(空虚な instantiation を Found と
      主張しない)に反する。

      **結論(訂正・正直な評価)**: このセッションで確立した一連の
      成果(Eisenstein拡大の全分岐・DVR・具体的な長さの下限)は、
      Theorem 1.2 の証明が要求する`Vₙ`塔側の**本物の、再利用可能な
      構成要素**として無駄になっていない——「典型例」の`V_{n+1}/V_n`
      をこのセッションの道具で**具体的に実装できる**ことを示した点は
      価値がある。しかし`δₙ→0`という**Theorem 1.2そのものの主張**を
      非空虚に完成させるには、`Wₙ`側(`Vₙ⊗_V W`の正規化、`W`の
      almost étale性を経由する)の平行した構成が要り、これは
      §2-4と**同根**(mathlibに無いHochschild cohomology・almost
      mathematics)の壁に行き着く——「`hrec`への最後の統合」は
      当初想定したような小さな仕上げ作業ではなく、**§2-4と同規模の
      新規理論構築を要する**ことが分かった。次回以降は、この
      `Vₙ`側の道具(`falt1_theorem12_length_differentIdeal_eq_
      length_quotient_n_add`系)を明示的に「Theorem 1.2の`V_n`塔
      構成部品」として文書化した上で、`Wₙ`側(almost étale)への
      本格着手が必要になる時点まで、この事実を正直に記録しておく。

      ★2026-09-05、続けて**§2への迂回路を検討した**: `[Falt1]
      Definition 2.1`(`isAlmostEtaleCovering`)は`Found/Falt1/
      AlmostEtale.lean`に既に`IsAlmostEtaleCovering`として sorry無く
      formalizeされている(commit `d1f92c36`、以前のセッション)が、
      **skeletonの`isAlmostEtaleCovering`は依然として`Interface.
      AlmostEtaleSetup.isAlmostEtale : Prop`という posit を指しており、
      本物の`IsAlmostEtaleCovering`に差し替えられていない**——差し替え
      て`AlmostEtaleSetup.example`の`isAlmostEtale := True`という
      **空虚な**instantiationを、本物の非空虚な witness に置き換え
      られれば、Definition 2.1(§2の4項目の1つ)が正直に Found になる
      可能性があると考え、その non-vacuous witness の構成に着手した
      (`Found/Falt1/AlmostEtale.lean`自身のdocstringが「A=Bの具体例
      でのwitness構成は未完成」と記録していた箇所)。

      **診断**: `A=B=R`(恒等拡大)の場合、`awayAlgebra`が構成する
      `Algebra (Localization.Away 1) (Localization.Away (algebraMap
      R R 1))` instanceと、mathlibが自動的に見つける標準の
      `Algebra.id`(自己代数)instanceが**衝突**する——`Localization.
      awayMap (algebraMap R R) 1 = RingHom.id _`(`IsLocalization.
      ringHom_ext`で証明可能、確認済み)+`(RingHom.id S).toAlgebra =
      Algebra.id S`(`rfl`)という**2つの橋渡し補題は確立できた**が、
      これを`Module.Free`等の**instance引数**として埋め込まれた
      ゴールへ実際に反映させる(`▸`/`rw`が instance引数には直接効かない
      ため、`convert`や明示的な `@`付き項に頼る必要がある)段階で、
      このセッションで初めて遭遇する種類の詰まりに当たった
      (`Module.Free.of_equiv`系の補題は`LinearEquiv`の型引数解決で
      別のエラーになった)。

      **回避策の発見**: `A=B`という**自己拡大**が衝突の根本原因だと
      気づいた——`B:=R×R`(`A≠B`の本物の非自明拡大、`p:=(1,1)`が
      unit)を使えば、`awayAlgebra`が構成する`Algebra (Localization.
      Away 1) (Localization.Away (1,1))`と**競合する標準instanceは
      そもそも存在しない**(`Localization.Away 1`と`Localization.
      Away (1,1)`は無関係な抽象型なので、`Algebra.id`のような自動
      発見される代替が無い)。この場合`Module.Free`等は`IsLocalization.
      atUnit`(`R ≃ₐ[R] Localization.Away p`、任意のunit`p`に一般化
      された`atOne`)経由で`R×R`が`R`上free・finite・étaleであること
      (mathlibに標準的にあるはず)へ帰着できる見込みだが、**trace
      写像の条件(ii)・idempotentの条件(iii)を含む完全なwitnessは
      今回は完成しなかった**——`A=B`の衝突を回避する方針は見えたが、
      残る作業量(4条件すべての具体的な証明)を考慮し、無理に押し切って
      壊れたものを残すより、この診断結果(2つの橋渡し補題+`A≠B`回避策)
      を正直に記録し、次回への持ち越しとした(commitは無し、Lean
      ファイルへの変更も無し)。

      ★★2026-09-05、続けて`B:=Fin 2 → R`(`R×R`と本質的に同じだが
      Pi型として扱うことでmathlibの一般論が直接使える)で、**`Etale
      R (Fin 2 → R)`が実際に成り立つことをscratch fileで確認した**
      (`lean/ABC3/Found/Falt1/AlmostEteleWitnessTest.lean`として
      一時的に作成・`lake build`成功・6秒・確認後削除、commitはして
      いない):
      - `Algebra.FormallyUnramified.pi_iff`(`Mathlib.RingTheory.
        Unramified.Pi`)+`Algebra.FormallySmooth.pi_iff`(`Mathlib.
        RingTheory.Smooth.Pi`、この2ファイルは`KaehlerAux.lean`の
        推移的importには含まれておらず`lean_check`のREPLでは見えない
        ——実際の`lake build`でのみ検証できた)により、`FormallyUnramified
        R (Fin 2→R)`・`FormallySmooth R (Fin 2→R)`はそれぞれ`∀i,
        FormallyUnramified/FormallySmooth R R`(自明、`self`インスタンス)
        に帰着する。
      - `Algebra.FinitePresentation.pi`(`Mathlib.RingTheory.
        Finiteness.FinitePresentationLocal`、`Mathlib.RingTheory.
        FinitePresentation`とは**別ファイル**、要注意)で
        `FinitePresentation R (Fin 2→R)`も同様に閉じる。
      - `Algebra.Etale.iff_formallyUnramified_and_smooth`でこれらを
        束ね、`Algebra.Etale R (Fin 2 → R)`が確立できた。
      **残る作業**: (a) この`Etale R (Fin 2→R)`(および同様に閉じる
      はずの`Module.Free`・`Module.Finite`)を、`awayAlgebra`が構成
      する`Algebra (Localization.Away 1) (Localization.Away
      (algebraMap R (Fin 2→R) 1))`という**具体的なinstance**の下へ、
      `IsLocalization.atUnit`(`R ≃ₐ[R] Localization.Away p`、任意の
      unit`p`版の`atOne`)経由で**移送する**作業(`Algebra.Etale.
      of_equiv`等はBASE環が同じ場合のみを扱うため、BASE環自体が
      `R→Localization.Away 1`と変わる今回は、両側を同時に移送する
      追加の橋渡しが要る)、(b) trace写像の条件(ii)、(c) idempotentの
      条件(iii)——これらは今回未着手。`Etale`性の核心部分(mathlibの
      一般論だけで閉じる)が確認できたことは、Definition 2.1の
      non-vacuous witnessが**原理的に到達可能**であることの具体的な
      裏付けになった。

      ★2026-09-05、続けて上記(a)(BASE環の移送)に着手し、正確な
      API シグネチャを`#check`で確定させたが、**完全には閉じなかった**
      ——次回のために詳細を記録する:
      ```
      Algebra.Etale.baseChange (R A B : Type*) [...] [Etale R A] :
        Etale B (B ⊗[R] A)
      IsLocalization.Away.tensorRightEquiv {R} (S) [...] (r:R) (A)
        [IsLocalization.Away r A] : A ⊗[R] S ≃ₐ[S] Localization.Away (algebraMap R S r)
      IsLocalization.Away.tensorEquiv {R} (S) [...] (r:R) (A)
        [IsLocalization.Away r A] : S ⊗[R] A ≃ₐ[S] Localization.Away (algebraMap R S r)
      ```
      `Algebra.Etale.baseChange R (Fin 2→R) (Localization.Away 1)`で
      `Etale (Localization.Away 1) (Localization.Away 1 ⊗[R] (Fin 2→R))`
      は直接得られる(BASE環が`Localization.Away 1`に**正しく**なる)。
      しかし、これを`Localization.Away(algebraMap R (Fin2→R) 1)`に
      接続する`tensorRightEquiv`(`S:=Fin2→R,r:=1,A:=Localization.Away
      1`で`Localization.Away1 ⊗[R] Fin2R ≃ₐ[Fin2R] Localization.Away
      (algebraMap R Fin2R 1)`)は**`≃ₐ[Fin2R]`(右因子側)であって
      `≃ₐ[Localization.Away 1]`ではない**——欲しい向き(`tensorEquiv`、
      `S:=Localization.Away1`側を左因子・base双方にする)は
      `[IsLocalization.Away r A]`(`A:=Fin2R`)を要求するが、これは
      一般に**成立しない**(`algebraMap R (Fin2→R)`は単射だが全射でない
      ため、`Fin2→R`はそもそも`R`の局所化ではない)。すなわち、
      `tensorRightEquiv`/`tensorEquiv`のどちらも「欲しい向き
      (`Localization.Away 1`が base **かつ** 左因子)」を直接与えない
      ——`AlgEquiv.ofRingEquiv`で`tensorRightEquiv`の**土台となる
      RingEquiv**を取り出し、それが`awayAlgebra`の`algebraMap`
      (`Localization.awayMap`)と可換であることを別途示す経路は
      残っているが、これは結局「`awayAlgebra`のinstanceと自然な
      instanceを橋渡しする」という**当初の困難と同型の作業**に
      帰着する——無理に押し切らず、この正確なAPI調査結果(シグネチャ
      3つ)を次回への足場として記録し、今回はここで打ち切った
      (scratch fileはlake buildで検証後削除、commitは無し)。

      ★★★★★★★★★★2026-09-05、**戦略を転換して(a)を完全に解決した**。
      `AlgEquiv`/テンソル積の橋渡し(`tensorRightEquiv`等)を諦め、
      **`RingHom.Etale`(bare ring homの性質、`Algebra`インスタンスを
      一切参照しないため diamond がそもそも起こらない)のレベルまで
      降りる**戦略に転換した:`Localization.awayMap (algebraMap R
      (Fin2→R)) 1`(`awayAlgebra 1`が使う環準同型そのもの)を、`p=1`が
      単元であることから得られる2つの全単射(`IsLocalization.atUnit`)
      `ιR : R≃ₐ[R]Localization.Away1`・`ιB : (Fin2→R)≃ₐ[Fin2→R]
      Localization.Away(algebraMap R(Fin2→R)1)`を使って`ιB∘algebraMap
      R(Fin2→R)∘ιR.symm`に分解する式`heqmap`を、`AlgEquiv`ではなく
      **`RingEquiv.ofBijective`+`IsLocalization.ringHom_ext`**(局所化
      からの環準同型は台の乗法的集合上での値だけで一意に決まる、という
      原理)で直接示した。これで`RingHom.Etale.of_bijective`(同型は
      étale)・`RingHom.etale_algebraMap`(`etale_fin2`——`Fin2→R`が`R`
      上étaleなことは`Algebra.FormallyUnramified.pi_iff`・`Algebra.
      FormallySmooth.pi_iff`で各成分に帰着、自明)・`RingHom.Etale.
      stableUnderComposition`(合成安定性)を貼り合わせるだけで
      `Algebra.Etale(Localization.Away1)(Localization.Away(algebraMap
      R(Fin2→R)1))`(`awayAlgebra 1`のもとで)が**sorry無く完成した**。
      `Module.Free`・`Module.Finite`も同じ`heqmap`から、`ιB`を
      半線形同値`(Fin2→R)≃ₛₗ[ιR.toRingHom]Localization.Away(...)`に
      仕立てて`Module.Free.of_equiv`、`ιR`・`ιB`からの環準同型の可換
      四角形を`he`として`heqmap`から直接示して`Module.Finite.
      of_equiv_equiv`(`≃ₗ`を経由しない代数レベルの移送、`Module.
      Finite.equiv`は同一環上の`≃ₗ`しか受け付けないため使えなかった)
      に渡すことで完成させた。3本ともまずscratch file(`AlmostEtele
      WitnessTest3.lean`、`lake build`で個別に検証)で確立してから、
      実プロジェクトの`Found/Falt1/AlmostEtale.lean`(既存の
      `awayAlgebra`をそのまま使うよう書き換え)へ移植し、**移植後も
      1回で`lake build`が通った**(`awayAlgebra`が`@[reducible]`な
      ため、scratch版の生の`(...).toAlgebra`表記との defeq が
      自動的に効いた)。プロジェクト全体の`lake build`(6590 jobs)も
      sorry無く成功、`node tools/check.mjs --brief`もNG13件(既存分)
      で不変を確認。`awayOne_fin2_etale`・`awayOne_fin2_freeFinite`
      として commit(`tools/lean-idioms.md` #44に手法を記録)。

      **意味**: `Definition 2.1`(`IsAlmostEtaleCovering`)の non-vacuous
      witness、条件(i)の3点(`Module.Free`・`Module.Finite`・`Algebra.
      Etale`)が`A:=R`・`B:=Fin2→R`・`p:=1`について**すべて完成した**。
      残るは条件(ii)(trace写像が`B`を`A`へ写す)・条件(iii)(idempotent
      `p^n e_{B/A}`が`B⊗_AB`の像に入る)の2点のみ——これが埋まれば
      §2(4項目)のうち`Definition 2.1`が初めて non-vacuous に Found と
      なる、§2-4(11項目)の総ブロック状態への最初の穴になる。

      ★★★★★★★★★★★2026-09-05(続)、**残る条件(ii)(trace)・条件(iii)
      (idempotent)も同じセッションで完成させ、`Definition 2.1`の
      non-vacuous witness を`awayOne_fin2_isAlmostEtaleCovering`として
      完全に組み立てた**。

      条件(ii)は**驚くほど自明だった**——`p=1`なので`algebraMap R
      (Localization.Away 1)`自体が全単射(条件(i)の`hbijR`と同じ事実)、
      すなわち**どんな元にも原像がある**。trace の値を計算する必要すら
      無く、単に`hbijR.2`(全射性)を trace の値に適用するだけ。

      条件(iii)は`p^n=1`(`p=1`なので)・`1•x=x`より`∃e, diagonalCompare
      1 e = elem Ap Bp`に帰着する。ここで鍵となる発見: **`elem Ap Bp`
      の値を計算・特定しようとせず、`diagonalCompare 1`の全射性だけを
      示せば良い**(`elem`は`Exists.choose`で非構成的に定義されており、
      具体的な計算式がmathlibに無いため、値を特定する経路は一意性補題
      を要し大掛かりになる——全射性経由はこれを完全に回避する)。
      `diagonalCompare p`の純テンソルでの値を`diagonalCompare_tmul`
      (`f0 := algebraMap B Bp`を両成分に施すだけ、と直接計算で確認)
      として切り出し、`f0`が全射なら`diagonalCompare p`も全射
      (`TensorProduct.induction_on`の3ケースで、`tmul`ケースだけ`f0`の
      全射性に帰着させる)ことを`diagonalCompare_surjective_of_
      algebraMap_surjective`として一般に証明した。`p=1`の場合`f0`は
      単元での局所化写像なので全単射(条件(i)の`hbijB`と同じ)——これを
      適用して`elem Ap Bp`自体への原像を直接引くだけで条件(iii)が
      閉じた。

      最後に`awayOne_fin2_isAlmostEtaleCovering : IsAlmostEtaleCovering
      (A:=R)(B:=Fin2→R)(1:R)`として(i)(ii)(iii)を`refine ⟨_,_,_,?_,?_⟩`
      で貼り合わせ、`R:=ℤ`での具体例(`example : IsAlmostEtaleCovering
      (A:=ℤ)(B:=Fin2→ℤ)(1:ℤ) := ...`)で非空虚性を確認した。まずmcp
      REPL(`lean_check`)で全パーツを個別に検証してから実ファイルへ
      移植し、**移植後も1回で`lake build`が通った**。プロジェクト全体
      の`lake build`(6590 jobs)・`node tools/check.mjs --brief`
      (NG13件、既存分で不変)も確認済み。手法を`tools/lean-idioms.md`
      #45に記録(「`Exists.choose`の元への到達は、値を計算せず経由する
      写像の全射性だけで示す」という一般に使える戦略)。

      **結論**: `Definition 2.1`(`IsAlmostEtaleCovering`)は、定義の
      formalization(sorry無し、既に完成済みだった)に加えて、**その
      定義を満たす具体的な非空虚な実例が実在することまで含めて完全に
      Found になった**——`Falt1 Chapter I`の13項目のうち、`Lemma 1.1`
      に続く**2件目**。§2-4(11項目)の総ブロック状態に初めて風穴が
      開いた。`Theorem 2.2`-`2.4`・§3・§4は引き続き Hochschild
      cohomology を要し未着手のまま。

      ★★★★★★★★★★★★2026-09-05(続々)、**同じ手法を`B:=Fin2→R`に限らず
      一般化し、以前「次のセッションへ持ち越す」としていた`A=B`
      (恒等拡大)の instance diamond も解消した**。`awayOne_fin2_*`の
      証明本体は`B`が`Fin2→R`である具体的な事実を一切使っておらず、
      `[Algebra.Etale A B][Module.Finite A B][Module.Free A B]`という
      **仮定**だけで同じ議論が通ることに気づいた——`awayOne_etale_of_
      etale`・`awayOne_freeFinite_of_etale`・`awayOne_trace_of_unit`・
      `awayOne_idempotent_of_etale`・`awayOne_isAlmostEtaleCovering_
      of_etale`として一般化し、`B:=A`(恒等拡大)を代入するだけで
      `IsAlmostEtaleCovering (A:=A)(B:=A)(1:A)`が得られることを確認
      した(`elem_self`の docstring が指していた本来の疑問への回答)。

      一般化の途中で新しい罠に当たった: `{A B : Type*}`と素朴に書くと
      `A`・`B`が**別々の**universe metavariableになり、`RingHom.Etale.
      stableUnderComposition`(`{R S T : Type u}`と**単一**のuniverseを
      要求)が意味不明な type mismatch で失敗する——`universe u`を
      明示的に宣言して`{A B : Type u}`と揃えることで解決した(`tools/
      lean-idioms.md` #46)。`Module.Free A B`/`Module.Finite A B`も
      `Algebra.Etale`単独からは従わない(開はめ込み等の反例がある)ため
      別途仮定として必要——`hsmul`の`show`が`rfl`で閉じなくなった箇所
      は`Algebra.smul_def`を`rw`で挟むことで解決した。

      移植後、実ファイルでも1回で`lake build`が通り、プロジェクト全体
      の`lake build`(6590 jobs)・`node tools/check.mjs --brief`
      (NG13件、既存分で不変)も確認済み。この一般化自体は新しい
      Found項目を生むものではない(`Definition 2.1`は既にFound)が、
      「`p=1`(単元)での almost étale covering は、étale・finite・free
      でありさえすれば任意の`B`について成立する」という、以前
      「解決困難」と記録していた具体的な gap を正式に閉じた、
      infrastructure としての確実な前進。

      ★★★★★★★★★★★★★2026-09-05(続々々)、**`p` を単元に限らない、
      真に非退化な non-vacuous witness を完成させた**——これは
      `Definition 2.1` の witness としてこれまでで最も意味のある形。

      鍵となった2つの発見: (1) **`Algebra.FormallyUnramified.elem`の
      一意性**(`elem_unique_of_props`)——`Exists.choose`で非構成的に
      定義された idempotent だが、定義性質(annihilate `1⊗s-s⊗1`・
      augment to `1`)を満たす元は実は一意であることを、`elem_absorb_
      of_prop`(定義性質(1)を満たす`t`について任意の`x`が`x*t=(μx⊗1)*t`
      という「吸収」性質を持つ)を経由して初めて mathlib に無かった形で
      証明した。(2) この一意性を武器に、**`diagonalCompare p (elem A
      B) = elem Ap Bp` を任意の`p`(単元である必要が無い)について
      証明した**(`diagonalCompare_elem_eq`)——`elem Ap Bp`の定義性質
      (1)を`Bp`の**全ての元**について確認する必要があり、`f0(B)`⊆
      `Bp`だけでは足りないが、`Z:={s'|(1⊗s'-s'⊗1)*t=0}`が**部分環**
      であること(`Zclosed_add`・`Zclosed_mul`)と、`Away`局所化の
      生成元`π:=f0(p)`の逆元も`Z`に入ること(`Zclosed_inv`、`(π⊗π)`で
      割ってから可逆性でキャンセルする論法)を示し、`IsLocalization.
      Away.surj`(`Bp`の任意の元は`f0(a)*π⁻ⁿ`の形)で全域に拡張した。

      この2つを使って、**「`B`が`A`上(古典的な意味で)étale・finite・
      freeでありさえすれば、任意の`p`についてalmost étale covering」
      という完全に一般の定理**(`isAlmostEtaleCovering_of_etale_
      general`)を組み立てた。条件(i)は`RingHom.Etale.propertyIsLocal.
      localizationAwayPreserves`(mathlibの「局所的な性質」framework、
      étale性の`Away`局所化保存性)・`Module.Finite.of_isLocalization`・
      `Module.free_of_isLocalizedModule`という、p=1witnessの`RingEquiv`
      経由の議論より遥かに直接的な道具で閉じた——`p`が単元かどうかに
      一切依存しない、はるかにクリーンな証明になった。条件(ii)は
      `Algebra.trace_localization`(mathlib既存)一発。条件(iii)は
      `p^n•`を`diagonalCompare`の`A`-線形性で外に出し、(2)を適用する
      だけ。**真の非単元素元`p:=5`・`B:=Fin2→ℤ`**での具体例
      (`5`は`ℤ`の単元ではない)で非空虚性を確認した——`p=1`退化ケース
      とは質的に異なる、genuinely non-degenerate な witness。

      移植後1回で`lake build`が通り、プロジェクト全体の`lake build`
      (6590 jobs)・`node tools/check.mjs --brief`(NG13件、既存分で
      不変)も確認済み。手法を`tools/lean-idioms.md` #48に記録。

      **意味**: `Definition 2.1`は今や「古典的な不分岐拡大は(任意の
      `p`について)almost étale coveringである」という、Faltings の
      理論の健全性を裏付ける非自明な事実まで込みでFoundになった。
      ただし`Theorem 1.2`本体(`W`が**分岐した**almost étale covering
      の場合の`δₙ→0`)には直結しない——不分岐`W`を選ぶと`Wₙ`も
      (étaleがbase changeで保たれるため)不分岐のままとなり
      `δₙ≡0`と自明になってしまう一方、Faltings の理論が本来扱いたい
      のは分岐したW(`Theorem 2.2`-`2.4`のHochschild cohomologyが
      要る場合)であるため。それでも、`elem`の一意性・自然性という
      道具自体は almost 数学の一般論として再利用可能な資産。
   β-(d+1)(δ_n-δ_{n+1})`、`β=min{1,δ_n/(d+1)}`)を整理した形
   `δ_{n+1}≤δ_n-min{1,δ_n/(d+1)}/(d+2)` から `δ_n→0` を、`V_n`・`W_n`
   の具体的構成に一切依存しない**純粋な実数列の不等式**として抽出・
   完全に証明した。2段の議論(δ_n≥d+1が続く限り固定量ずつ減るので
   有限回で脱出/以降は比 `1-1/((d+1)(d+2))∈(0,1)` の等比減衰)を
   `squeeze_zero`+`tendsto_pow_atTop_nhds_zero_of_lt_one` で結んだ。
   ★残るのは、この不等式の**前提**(`hrec`)を実際の `V_n`・`W_n` の
   `Module.length` から導く部分(3bの一部・元の「劣加水性・完全列」
   評価)——これは 3b と一体で、`V_n` 族の具体形(3c)が定まってから
   でないと着手しづらい。

**§2-4(11項目)**: ★★★★★2026-09-04、`Definition 2.1`(almost étale
covering の定義)は `Found/Falt1/AlmostEtale.lean` で sorry 無く
formalize できた(`IsAlmostEtaleCovering`、commit `d1f92c36`)——
「projective・trace・idempotent」は個別には mathlib の道具(dual basis
の補題・`Algebra.trace`・`Algebra.FormallyUnramified.elem`)で書け、
特に「`B⊗_AB` の像」を捉える比較写像は局所化とテンソル積の可換性の
同一視を経由せず `Algebra.TensorProduct.lift` で直接構成できることを
発見した(想定より遥かに単純だった)。

★★★★★★★続けて `Theorem 2.2`(nilpotent ideal 上の lifting)の証明全文
(物理p.7=印字p.260)を260dpiで精読した。証明は: (1) `B` の almost
射影性から任意の `ε>0` で `p^εφ` を `A` 加群写像として持ち上げる、
(2) 積の非可換性の障害
`b_0(p^εφ_ε(b_1b_2)-φ_ε(b_1)φ_ε(b_2))b_3` が **Hochschild
cohomology** `H^2(B/A,I)` の類を定め、これが `m` で零化される
(remark 2.1(v)、証明は「同様の議論」とだけ書かれ本文中では未証明)
ことから `ε` を倍にして消せる、(3) 持ち上げの一意性は `H^1(B/A,I)`
(これも `m` 零化)まで、`C` が `p` 捩れ無しなのでちょうど一意に決まる、
(4) `mB→C` への拡張は `p^ε` の2乗関係から。

**結論**: Definition 2.1 とは質的に異なり、`Theorem 2.2`-`2.4`・`§3`・
`§4` は **Hochschild cohomology(`B⊗_AB⊗_A⋯⊗_AB→M` の複体としてremark
(v)で定義される)の理論そのもの・その almost 消滅性・一般の変形理論
(`H^2`=障害・`H^1`=分類)** という、mathlib に存在しない独立した理論を
要する。Definition 2.1 で見つかった「見た目より簡単だった」という
僥倖はここでは再現しない見込み——引き続き完全に未着手。

★★★★★★★★★2026-09-05、**`Theorem 2.4`(物理p.7-8=印字p.260-261)を
260dpiで精読し、(i)(ii)で難度が質的に違うことを発見した**——次回への
具体的な足場として記録する。

- **(i)**(`Ω_A⊗_AB→Ω_B` が almost isomorphism)の証明は「証明はいつも
  通り Hochschild cohomology `H^*(B/A,M)` を使う」の一文で済まされて
  おり(`Theorem 2.2`と同じ「同様の議論」型の省略)、これは Definition
  2.1(v)(remark、H^*(B/A,M)の定義そのもの)を前提にした一般論——
  引き続きブロックされたまま。
- **(ii)**(有限群`G`の semilinear 作用、`m` が `H^i(G,M)`(`i>0`)・
  `M^G/tr_G(M)` を零化する)の証明は、**驚くほど自己完結的で
  Hochschild cohomology を要しない**——「`m` が `A/tr_{B/A}(B)` を
  零化することを示せば十分」に帰着し、これは Definition 2.1 直後の
  remark (iii)(物理p.6末尾、印字p.259)そのもの:「`p^ε・e_{B/A}=
  Σxᵢ⊗yᵢ`(条件(iii)の witness)なら、**任意の`b∈B`について
  `p^ε・b=Σtr_{B/A}(b・xᵢ)・yᵢ`**」という恒等式一発で閉じる
  (`b:=1`とすれば`p^ε=Σtr(xᵢ)yᵢ∈tr_{B/A}(B)・B`)。この remark
  (iii)自体は「`e_{B/A}`が(honestに étale な `Bp/Ap` 上で)trace形式
  の双対基底を与える」という古典的な事実の系——**mathlibには`elem`と
  trace形式の双対基底性を結ぶ補題が無い**(`traceForm_nondegenerate`
  等は体の拡大専用、一般の有限étale環拡大向けが無い)ため、この
  remark (iii) 自体を独立に証明する必要がある(未着手・次回の具体的
  な最有力候補)。

  **ただし** `Theorem 2.4(ii)`の**主張そのもの**(`H^i(G,M)`の消滅)を
  完成させるには、remark (iii) に加えて**群コホモロジー**
  (mathlibの`groupCohomology`名前空間、restriction-corestriction /
  transfer 論法を`|G|`の代わりに`tr_{B/A}`の almost 全射性で置き換える
  一般化)が別途要る——これは remark (iii) 単体より一段大きい、しかし
  Hochschild cohomology(未証明のH^2消滅を要する`Theorem 2.2`-`2.3`・
  `(i)`)ほど絶望的ではない**中間的な難度**の項目として記録する。

  **意味**: `Theorem 2.4`は`(i)`と`(ii)`で難度が質的に異なる——`(ii)`
  の remark (iii) 部分(自己完結的な trace 恒等式)は、`Definition
  2.1`の witness で見つかった僥倖と同種の「見た目より簡単」な部分を
  含む可能性がある。§2-4(11項目)全体が一様にHochschild cohomology
  でブロックされているという以前の評価は、`(ii)`単体について**部分的
  に訂正**する——ただし完全な`Theorem 2.4`(群コホモロジーの完成まで
  込み)を Found にするには、なお相応の追加作業(群コホモロジー
  infrastructure の構築)が要る。

  ★★★★★2026-09-05(続)、**remark (iii) の恒等式(`p^ε・b=Σtr(b・xᵢ)yᵢ`)
  を、`elem`の一意性で使った「annihilation性質だけから出発する」流儀で
  部分的に証明した**。`Tr1map R S : S⊗_RS→ₗ[R]S`(`x⊗y↦tr(x)•y`)を
  定義し、`Ψ_s(t):=Tr1map((s⊗1)*t)`が`s・Tr1map(t)`に等しいことを
  annihilation性質(swap: `(1⊗s)*t=(s⊗1)*t`)だけから示した
  (`Ψ_eq_smul_Tr1`)。これで remark (iii) の恒等式は**`Tr1map(elem
  R S)=1`という1個のスカラー等式**に完全に帰着した。

  この最後のスカラー等式(「trace form の非退化性」の系、分離拡大論の
  古典的事実)は、Fin2→Rの具体例・R自身(恒等拡大)の両方で手計算では
  `1`になることを確認したが(rankではなく常に`1`)、一般の証明は
  mathlibに欠けている(`traceForm_nondegenerate`系は体の拡大専用)。

  ★★★★★★★★★★2026-09-05(続々)、**この最後のスカラー等式
  `Tr1map(elem R S)=1`を、基底を用いた添字計算で完全に証明した**。
  `S`(`R`上有限自由)の基底`{bᵢ}`で`elem R S=Σᵢbᵢ⊗aᵢ`と分解する
  (`tensorDecompEquiv`、`TensorProduct.congr`+`TensorProduct.
  finsuppScalarLeft`で`S⊗_RS≃(ι→₀S)`を構成)と、annihilation性質から
  構造定数`b.repr(bⱼbᵢ)k`を経由した関係式`bⱼaₖ=Σᵢ(b.repr(bⱼbᵢ)k)aᵢ`が
  出る(`tensorDecompEquiv_extract`)。これと augmentation(`Σbᵢaᵢ=1`)・
  trace の行列表示(`Algebra.trace_eq_matrix_trace`+`Algebra.
  leftMulMatrix_eq_repr_mul`)を、`S`の可換性(構造定数の対称性
  `b.repr(bⱼbᵢ)k=b.repr(bᵢbⱼ)k`、`bⱼbᵢ=bᵢbⱼ`から)で束ねる
  (`Finset.sum_comm`)と、**`Tr1map(elem)`と`1`が同じ二重和の並べ替え**
  であることが分かり、証明が閉じた(`Tr1map_elem_eq_one`)。

  さらに、`Tr1map`と`diagonalCompare`の可換性(`lmul'_diagonalCompare`
  と同型の議論、`Algebra.trace_localization`経由、`Tr1map_
  diagonalCompare`)を組み合わせ、**`Theorem 2.4(ii)`の remark(iii)
  恒等式そのものを完全に証明した**(`remark_iii_trace_identity`)——
  `IsAlmostEtaleCovering`の条件(iii)の witness `e`から、Faltings が
  主張する`p^n・b=Σtr_{B/A}(b・xᵢ)yᵢ`を導く。`algebraMap B Bp`の単射性
  (Faltings 自身が「全ての環で`p`による乗法は単射」と明記している
  標準仮定)を経由する。

  移植後、実ファイルで1回で`lake build`が通り(警告2件のみ、リント
  修正で解消)、プロジェクト全体の`lake build`(6590 jobs)・`node
  tools/check.mjs --brief`(NG13件、既存分で不変)も確認済み。

  **意味**: `Theorem 2.4(ii)`の証明の核心(remark(iii)恒等式)が
  sorry無く完成した。残るは主張本体(`H^i(G,M)`の消滅、群コホモロジー
  ——mathlibの`groupCohomology`名前空間+transfer論法の一般化)のみ
  ——これは`Theorem 2.2`-`2.3`(未証明のHochschild cohomology H^2消滅
  を要する)とは質的に異なり、既存のmathlibインフラ(群コホモロジーは
  既に存在する)を組み合わせる作業に近い、より見通しの良い次回の
  課題として記録する。

  ★2026-09-05(続々々)、**`groupCohomology`名前空間を実際に調査した**
  (正しい import は`Mathlib.RepresentationTheory.Homological.
  GroupCohomology.Basic`——`Homological`が経路に要る、当初のガイダンス
  の経路は誤り)。`groupCohomology A n : ModuleCat k`(`A:Rep k G`、
  `k`線形な`G`表現のコホモロジー、`HomologicalComplex.homology`で
  定義済み)・`groupCohomology.H0/H1/H2`(低次数の具体形)は確認できた
  が、**「`|G|`(あるいは一般に transfer/corestriction 論法での
  重み)が`H^i(G,M)`(`i>0`)を零化する」という標準的な事実
  (restriction-corestriction合成が`|G|`倍になる、という古典的な
  transfer 論法)は mathlib に見当たらなかった**(`groupCohomology`
  名前空間 241 件を機械的に確認、それらしい補題名は無し)。この
  transfer 論法自体を`inhomogeneousCochains`(コチェイン複体、
  `ModuleCat`圏論的な取り扱い)のレベルから構築する必要があり、
  これは homological algebra の新しい API 領域(`ModuleCat`・
  `HomologicalComplex`の圏論的取り扱い)を要する、独立した規模の
  作業と判断した——`remark_iii_trace_identity`(今回完成)を
  `tr_{B/A}`版の transfer 写像の「重み」として使う、という数学的な
  筋道は明確になったが、それを実装する homological algebra の配管
  自体は次回への持ち越しとする。

  ★★★★★★★2026-09-05(続々々々)、**`Theorem 2.4(ii)`の証明本文が
  要求する「ノルムを両辺に適用する」ステップを、mathlibの`Ideal.
  relNorm`(相対イデアルノルム、`Mathlib.RingTheory.Ideal.Norm.
  RelNorm`、`IsDedekindDomain`+`Module.IsTorsionFree`な環同士)で
  実際に完成させた**(`trace_ideal_pow_mem_traceIdeal`)。
  `remark_iii_trace_identity`(`b:=1`)で得た`B`側のイデアル関係
  `algebraMap A B(p^n)∈tr_{B/A}(B)・B`を`Ideal.relNorm_algebraMap`
  (`relNorm(I.map(algebraMap)) = I^finrank`)で両辺に適用し`A`側へ
  引き戻すと、`p^{n・finrank}∈tr_{B/A}(B)^{finrank}⊆tr_{B/A}(B)`
  (`Ideal.pow_le_self`)が出る——これは Faltings の証明の最後の一文
  「`N_{B/A}`を両辺に適用して結論を導く」の**実際の数学的内容**に
  正確に対応する。`n`は条件(iii)で任意に選べるため、`n・finrank`
  (`finrank`個おきの指数)で`m`零化を実現する——`finrank=1`でない
  限り厳密には「全ての`n`」ではなく「`finrank`の倍数の`n`」だが、
  イデアルの吸収性(`p^{finrank}∈a`⟹任意の`m≥finrank`について
  `p^m∈a`)により実用上は十分な強さ。

  移植後、実ファイルで1回で`lake build`が通り(新規import`Mathlib.
  RingTheory.Ideal.Norm.RelNorm`が`IsDedekindDomain`系の依存関係を
  引くため初回コンパイルはやや長め)、プロジェクト全体の`lake build`
  (6590 jobs)・`node tools/check.mjs --brief`(NG13件、既存分で不変)
  も確認済み。

  **意味**: `Theorem 2.4(ii)`の証明本文が明示的に要求する3つのステップ
  (remark(iii)の trace 恒等式・ノルムの適用・`A/tr_{B/A}(B)`のalmost
  消滅)のうち、**最初の2つ(remark(iii)恒等式・ノルム適用)が両方とも
  sorry無く完成した**。残るは(a)`Ideal.pow_le_self`の`finrank`ギャップ
  を厳密に埋める(または「finrank個おき」で十分だと割り切る)、
  (b)`Theorem 2.4(ii)`が実際に主張する`H^i(G,M)`(`i>0`)への一般化
  (群コホモロジーのtransfer論法、`groupCohomology`名前空間の欠落分の
  構築)——この(b)が依然として最大の残工程。

★★★★★★★★★★★★★★★2026-09-05(続々々々々)、**重大な訂正**: これまで
複数回、「remark 2.1(v)(`Theorem 2.2`-`2.4`が使う「`m`がHochschild
cohomologyを零化する」という事実)はFaltings自身が本文中で証明せず
外部参照に頼っている」と報告してきたが、**これは誤りだった**。原文
(物理p.6=印字p.259)を260dpiで再確認したところ、remark (v) は
「`e_{B/A}=Σxᵢ⊗yᵢ`が`B⊗AB`の元だったなら、`b₀⊗b₁⊗⋯⊗b_{n+1}↦
Σxᵢy_ib₀⊗b₁⊗⋯⊗b_{n+1}`という null-homotopy が得られ、Hochschild
cohomology は消滅する。The same argument gives that for B almost
etale over A m annihilates the Hochschild cohomology in positive
degrees.」という**完全な、標準的な**構成を与えている——分離代数論の
古典的事実(バー分解の縮約ホモトピー)そのもので、Faltings は何も
省略していない。

この事実(の honest な場合、`p`が単元である必要が無い一般の formally
unramified拡大)を、mathlibの`CategoryTheory.Abelian.Ext`(導来圏経由
の一般Ext理論、`Mathlib.Algebra.Homology.DerivedCategory.Ext.*`)を
使って**証明した**(`hochschild_ext_eq_zero`)。`elem R S`の
annihilation性質から`S`が`S⊗_RS`-加群として`μ:S⊗_RS→S`の**切断**
(`hochSection`、`s(b):=(b⊗1)*elem`)を持つことを示し(`S⊗_RS`線形性
の検証が核心、swap annihilation性質`one_tmul_mul_elem`を使う)、これは
「`S`が`S⊗_RS`-加群として射影的」を意味する(`hochModule_projective`、
`Module.Projective.of_split`経由)。射影加群からの`Ext`は正の次数で
消えるという一般論(`Ext.eq_zero_of_projective`、mathlib既存)一発で、
`HH^n(S/R,M):=Ext^n_{S⊗_RS}(S,M)`が`n>0`で消えることが出た。

lake build(プロジェクト全体、6590 jobs、新規importが
`DerivedCategory.Ext`系の大きな依存関係を引くため初回はやや時間が
かかる)・`node tools/check.mjs --brief`(NG13件、既存分で不変)で
確認済み。

**意味**: remark (v) が Faltings の言う通り**成立する完全な定理**で
あることを実際に証明で実証した——以前の「未証明」という評価は
正式に撤回する。ただし今回証明したのは remark (v) の **honest な
場合**(`S`が`R`上honestに formally unramified、`elem`のannihilation
性質がexactに成り立つ場合)のみ——`Theorem 2.2`-`2.4`が要求する
**almost**な場合(`B`が単に almost étale、`p^ε elem`のみ`B⊗AB`に
ある場合)への一般化は、`remark_iii_trace_identity`と同型の「局所化
を経由してinjectivityで降ろす」議論を要し、`hochSection`が almost
idempotentでは`S⊗_RS`線形にならない(annihilation性質がexactに
成り立たないため)という技術的な障壁を今回確認した——次回への
具体的な課題として記録する。

★★★★★2026-09-05(続々々々々々)、**上記の技術的障壁を解消した**:
`diagonalCompare p : B⊗_AB → Bp⊗_ApBp`(`Ap:=Localization.Away p`・
`Bp:=Localization.Away(algebraMap A B p)`)が、`B`が`A`上free・
`algebraMap B Bp`が単射(Faltings自身が終始置く「p-torsion-free」
標準仮定)という2条件だけから**単射**であることを証明した
(`diagonalCompare_injective`)——`Module.Flat`の`rTensor`/`lTensor`
`preserves_injective_linearMap`を2回使い`φ:=(rTensor Bp f₀)∘(lTensor
B f₀)`の単射性を得て、`IsLocalization.moduleTensorEquiv`(「両成分が
既に局所化された`Ap⊗_ApBpBp`から、さらなる局所化`Ap`を消して
`A⊗Bp⊗Bp`に戻す」という mathlib の同型)経由で`diagonalCompare`を
`e.symm∘φ`と分解するのが鍵。

これにより、`Ap,Bp`レベルでのhonestなannihilation(`elem`の性質)を
`diagonalCompare`で引き戻し、単射性で`B⊗AB`へそのまま降ろせる
(`almost_swap_annihilate`・`almost_swap_mul_eq`)。augmentation側も
同様に`algebraMap B Bp`の単射性で`μ(w)=p^n`(`B`の中で厳密に)まで
降ろせた(`almost_swap_augment`)。`hochSection`自体を任意のwitness
`w`とそのswap性質だけから構成できるよう一般化し(`hochSectionOfWitness`)、
上記を組み合わせて`hochSectionAlmost_augment`——「`S`が`S⊗AS`の
"almost direct summand"(`μ∘s=`「`p^n`との掛け算」)」という remark(v)
のalmost版主張を、`IsAlmostEtaleCovering`のみ(`B/A`自体がhonestly
unramifiedである必要なし、真に一般の almost 設定)から**完全に
形式化**した。lake build(6590 jobs)・`node tools/check.mjs --brief`
(NG13、既存と不変)で確認済み(コミット参照)。

★★★★★★2026-09-05(続々々々々々々)、**残っていた最後の1段も閉じた
——remark 2.1(v)の almost 版が完成した**(`hochschild_ext_almost_zero`)。

「almost projective ⟹ almost Ext消滅」の論法(`s:B→T`・`μ:T→B`
(`T:=B⊗AB`)について`s≫μ=τ•𝟙`なら、`Ext^{k+1}(T,M)=0`(`T`は`T`上
射影的)を経由して`τ`が`Ext^{k+1}(B,M)`を零化する)は、mathlib の
`Ext.mk₀`(射から`Ext ⋯ 0`)・`Ext.mk₀_comp_mk₀_assoc`(合成の関手性)
・`Ext.smul_comp`/`Ext.mk₀_smul`(`R`-線形圏での`Ext`の`R`-加群構造)
・`Ext.eq_zero_of_projective`だけで**そのまま書ける**ことが分かった
(`ext_smul_eq_zero_of_almost_split`、圏論の言葉だけの一般補題として
切り出した)。

**配管上の落とし穴を1つ発見**: `Ext X Y n`の`R`-加群構造は
`Mathlib/Algebra/Homology/DerivedCategory/Ext/Linear.lean`にあり、
`Ext/EnoughProjectives.lean`(honest 版で使っていた import)からは
**推移的に import されない**。これに気づかない間は
`Module T (Ext S M (k+1))`のインスタンス探索が静かに失敗し続け、
しかも`#check @CategoryTheory.Abelian.Ext.<未知の名前>`が
「unknown identifier」ではなく`HasSmallLocalizedHom`の探索
タイムアウトを返す(`Ext`が`def`なので generalized field notation
として解釈されるため)ので、**名前の存在確認自体が誤誘導される**。
`tools/lean-idioms.md #51`に登記した。

非空虚性の対照も付けた——真の非単元`p:=5`・`A:=ℤ`・`B:=Fin 2 → ℤ`で
仮定4つ(almost étale covering・`Module.Free`・`algebraMap B B[1/5]`の
単射性・条件(iii)の witness)がすべて実際に成り立つ。

**残るのは`Theorem 2.2`-`2.4(i)`の"定理本体"**——「Hochschild
cohomology が`m`で零化される」という**入力**は完成したので、次は
それを使う議論(nilpotent ideal に沿った lifting の障害類を実際に
`Ext²`の元として取り出す部分)の設計になる。

★★★**結論(正直な評価)**: `/goal` の 13/13 は、このセッションでの
継続作業だけでは現実的な時間内に到達できない規模の作業(Theorem 1.2
1件でさらに大きな作業、§2-4 は新規ライブラリ構築に近い規模)である
可能性が高い。ただし Lemma 1.1 で確立した一般補題(especially
`subsingleton_H1Cotangent_self`・`polynomialKaehlerSplit`)は Theorem
1.2 の帰納法でも再利用できる見込みがあり、**無駄になっていない**。

★★★★2026-09-04、**§3・§4 の主張そのものを初めて実際に読んだ**
(`lean/ABC3/Skeleton/Falt1/Section3.lean`・`Section4.lean`、既存の
Interface posit を精読、新規の 260dpi 目視は不要——既にファイル内に
逐語が写してあった)。これまで「§2 が Hochschild cohomology を要する
ので恐らく §3・§4 も同様に重い」という**推測**だったが、今回**§3・§4
自身の主張を具体的に確認**して評価を裏付けた:

- **§3(Theorem 3.1・3.2、2項目)**: `S⊗_R R_∞` の正規化が `R_∞` の
  almost étale covering であること——`Interface/Falt1/GoodReduction.lean`
  の `GoodReductionSetup` が既に `AlmostEtaleSetup.isAlmostEtale`
  (旧い posit 述語、**まだ `Found/Falt1/AlmostEtale.lean` の実際の
  `IsAlmostEtaleCovering` に差し替えられていない**)を参照する形で
  スケルトン化済み。証明そのものは Theorem 2.2-2.4(Hochschild
  cohomology による持ち上げ)の結果を使う——**§2 の壁がそのまま
  継承される**ことを確認した。
- **§4(Theorem 4.1-4.3・4.5、5項目)**: `H^i(Δ,R̂)` という Galois
  コホモロジー(`Δ` = Galois 群、`R̂` = 完備化された環)の計算、
  スペクトル系列 `E_2^{a,b}=H^a(Δ,H^b(I,R̂))⇒H^{a+b}(Δ,R̂)` の
  almost 退化、cup 積による Hodge-Tate 型の同型——これは**§2 より
  さらに深い**、p進 Hodge 理論の核心部分そのもの(周期環の Galois
  コホモロジー・スペクトル系列)。mathlib には(1)almost mathematics
  だけでなく(2)この種の Galois コホモロジー計算の対象となる周期環
  `R̂` の理論も存在しない——**2重に新規ライブラリ構築が要る**。

★結論(確定・更新): §2-4 のブロッカーは**単一の Hochschild
cohomology だけではなく**、§4 では独立にもう1段深い(周期環の
Galois コホモロジー・スペクトル系列)理論を要することが**具体的に
確認**できた——「推測」から「該当ページを読んで確認した事実」への
格上げ。これで §2-4(11項目)全てについて、着手前の偵察は完了した
と見なしてよい——次に必要なのは実際の mathlib 開発であり、この
セッション規模の作業ではない。

## 1. 構造(2026-09-04 確定——旧版の3つの未確認点をすべて解消した)

★★**Chapter見出しは"CHAPTER"の語を伴わない裸のローマ数字**("I""II""III"
が単独行に出るだけ)。原文 物理p.4(印字p.256)の(e)段落「... In §I we study
the case of good reduction ... In §II we construct the theory 𝒳*(X) ...
Finally in §III we treat the case of bad reduction ...」で著者自身が
3章構成を明言しており、実測(下表)と完全に一致する。

★★**項目の番号付けは "N.M. Kind" の順**(番号が先、種別が後)。LocProP等の
"Kind N.M"(種別が先)と**逆**。`tools/paper-items.mjs`はそのまま使えない
(専用の抽出スクリプトが必要——再現用ロジックは本ファイル末尾のメモ参照)。

★★**pageOffsetを訂正した(254→253)**: 物理p.1はJSTORの表紙(著者・出典
スタンプのみ、本文ページ番号なし)で、物理p.2が印字p.255(論文本文の
実際の1ページ目)。150dpiで物理p.1-3を目視し、物理p.3の見出し
"256 GERD FALTINGS"で確定した。★[LocProP]形式化時に登記した
`theorem_I_4_4.src`の`pdfPage := 17`は260dpi目視による直接確認だった
ため訂正の影響を受けない(訂正後の式 印字270-253=17 で無矛盾も確認済み)。

Chapter III §4 は末尾を確認した結果、**4.1 Theoremの1件のみ**(物理p.45、
印字p.298)で、直後に p.299 の BIBLIOGRAPHY が続き論文が終わる。
→ 総項目数は **30 で確定**(旧版の「30+」の「+」は不要)。

## 2. 全項目表(物理ページは pageOffset=253 で算出・機械抽出。個別の
   260dpi目視は §4.4 以外まだ)

### Chapter I(good reduction、4節・13項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し(OCR、地の文の空白喪失あり) |
|---|---|---|---|---|
| §1 Ramification theory for discrete valuation rings | 1.1 | Lemma | 4 | For any extension V⊂W... the natural map... |
| | 1.2 | Theorem | 5 | If V⊂W is any extension and 𝒥'V denotes the normalization... |
| §2 Almost unramified extensions | 2.1 | Definition | 6 | Suppose A is a ring, B an A-algebra. B is called an almost... |
| | 2.2 | Theorem | 7 | Suppose B=A+mB is an almost etale covering of A, C an... |
| | 2.3 | Theorem | 7 | Suppose I⊂A is a nilpotent ideal, B an almost etale covering... |
| | 2.4 | Theorem | 8 | Suppose B is an almost etale covering of A. |
| §3 Good reduction | 3.1 | Theorem | 10 | Suppose S is a normal finite R-algebra such that S[1/p]... |
| | 3.2 | Theorem | 12 | Suppose S is a finite normal torsionfree R-algebra such that... |
| §4 Differentials and cohomology | 4.1 | Theorem | 13 | (i) The map Ω_{R/V}⊗R̄→... induces almost isomorphisms |
| | 4.2 | Theorem | 15 | (i) H^i(Δ,R̂)≅Λ^i((R⊗_V V̄)^(-1)d)⊕(rest)... |
| | 4.3 | Theorem | 17 | There exists a Γ-equivariant functorial extension |
| | **4.4** | **Theorem** | **17** | **(i) H^i(Δ,R̂)≅Λ^i((R⊗_V V̄)^(-1)d)⊕(rest)...**★先行あり(縮小形) |
| | 4.5 | Theorem | 19 | (i) The spectral sequence |

### Chapter II(𝒳*(X) の構成、3節・8項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し |
|---|---|---|---|---|
| §1 Construction of 𝒢* | 1.1 | Theorem | 23 | (i) There exists a spectral sequence |
| | 1.2 | Theorem | 25 | (i) There exist spectral sequences |
| §2 The isomorphism with étale cohomology | 2.1 | Lemma | 27 | Any x∈X is contained in an open U⊂X which is a K(π,1). |
| | 2.2 | Corollary | 27 | For X as above, x∈X, Spec(𝒪_{X,x̄}⊗V K̄) is a K(π,1). |
| | 2.3 | Lemma | 28 | Suppose X is smooth over V, D⊂X a divisor with normal... |
| | 2.4 | Theorem | 30 | The transformation (X proper and smooth, D⊂X a divisor... |
| §3 Relations to Hodge cohomology | 3.1 | Theorem | 35 | Suppose X is proper and smooth over V, D⊂X a divisor with... |
| | 3.2 | Theorem | 36 | Suppose X is proper and smooth over V=W(k)(Witt vectors)... |

### Chapter III(bad reduction、4節・9項目)

| 節 | 項目 | 種別 | 物理p. | 内容の頭出し |
|---|---|---|---|---|
| §1 Commutative algebra | 1.1 | Lemma | 37 | There exists an e independent of n such that the normalization... |
| | 1.2 | Proposition | 38 | Suppose R has one system of good units for some (infinite)... |
| | 1.3 | Theorem | 39 | Suppose that R is a V-algebra of finite type, geometrically ir-... |
| §2 Global methods | 2.1 | Definition | 40 | A stable punctured curve of type (g,r) over an algebraically... |
| | 2.2 | Lemma | 41 | If V is a complete discrete valuation ring with fraction field K,... |
| §3 Rigid coverings | 3.1 | Definition | 41 | A V-map f:X→Y is called a rigid étale covering if |
| | 3.2 | Theorem | 42 | Suppose f:U.→Y is a rigid étale hypercovering, a coherent... |
| | 3.3 | Theorem | 43 | Suppose R is a V-algebra of finite type, smooth in characteristic... |
| §4 Intermediate cohomology | 4.1 | Theorem | 45 | Suppose X is a smooth proper K-scheme, D⊂X a divisor with... |

## 3. [LocProP]が直接引用する箇所(確認済み・実装済み)

| LocProP側 | Falt1側 | 状態 |
|---|---|---|
| [LocProP] Lemma 2.1 | Chapter I, Theorem 4.4, (i)(iv) | ★★物理p.17(印字p.270)で260dpi目視確認済み。実装済み(縮小形) |
| [LocProP] Lemma 2.3 | Chapter I(Leray-Serre + 前2つ、Theorem 4.3?) | 未確認・未着手 |
| [LocProP] 全体(§2の土台) | almost étale extension理論(§I-2) | 未確認・未着手 |

## 4. 着手順序(葉から)

leaf-first の原則により、他の項目に依存しない・短い証明から着手する。
候補(未検証・着手時に依存グラフで裏取りすること):
1. **Chapter I §1**(1.1 Lemma, 1.2 Theorem)—— almost étale extension理論
   の手前の純代数的補題で、他の Falt1 内項目への依存が最も少なそう。
2. **Chapter III §2.1**(Definition, stable punctured curve)—— 定義のみ
   なら証明義務が無く着手しやすい。
3. Chapter I §4.4 の縮小形を(i)-(vii)全体に拡張(LocProPの土台なので
   優先度は高いが、(ii)(iii)(v)(vi)(vii)は almost étale extension の
   深い性質を要する可能性があり、葉ではない懸念がある——着手前に
   `node tools/graph.mjs`で確認)。

## 5. 番号付け規則の再現手順(scratchpad/falt1-items.mjsの後継)

`scratchpad/falt1-full-index.mjs`(本セッションで作成、再現可能):
- 章境界: 行が正規表現 `^\s*(I|II|III)\s*$` に一致する単独行を探す
  (`CHAPTER`という語は原文に無い)。
- 節見出し: `^\s*(\d+)\.\s*[A-Z][a-z]` (空白喪失により単語がくっつくが
  見出し行自体は検出できる)。
- 項目: `^\s*(\d+)\.(\d+)\.\s*(Theorem|Proposition|Lemma|Corollary|
  Definition|Remark)`("N.M. Kind"の順、LocProPと逆)。
- 印字ページ: 走り込み見出し行 `NNN GERD FALTINGS`(偶数ページ)・
  `p-ADIC HODGE THEORY NNN`(奇数ページ)を全行から収集し、各項目の
  直前・直後で挟んで判定。
- 物理ページ = 印字ページ - 253(§1で確定した値)。

## 6. 未検証事項(次の個別着手時につぶすこと)

- 上表の物理ページ番号は§4.4以外、機械抽出のみ(260dpi目視未実施)。
  逐語照合(`原文 (Falt1 p.N):`)を使う前には必ず260dpiで該当ページを
  確認すること(JSTORスキャンは地の文の単語間スペースを失うが、定理の
  数式・番号は無事——Theorem 4.4の実測で確認済みの傾向)。
- 各項目が Falt1 内の他項目にどう依存するか(依存グラフの辺)は未調査。
  着手順序の葉判定は個別に`node tools/graph.mjs --owner Falt1`で
  裏取りすること。

関連: [[measure-mathlib-before-skeleton]] / `ResearchPaper/mathlib-gap.json`の
`locprop-perfectoid-hodge-tate`(この論文が主典拠)
