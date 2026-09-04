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

3c. **「非常に分岐した」`V_n` の族**そのものの形式化(具体例: `p^n`
    乗根と `1` のべき根を添加する塔、上の「典型例」段落)——まだ
    手つかず。抽象的な族の公理化(`Ω_{V_n/V_{n-1}}` が
    `(V_n/pV_n)^{d+1}` を商に持つ、という条件だけを構造として
    持たせる)なら 3a・3b より軽い可能性がある——ただし 3b が独立の
    length 計算を要すると分かった今、3c だけを先に進めても
    Theorem 1.2 全体の完成には直結しない(3b を経由しない限り
    `delta_tendsto_zero` の `hrec` に接続できない)。
4. ★★★★★2026-09-04、**完成した**(`delta_tendsto_zero`、commit
   `a9faa64e`)。長さの漸化不等式(上の逐語引用の通り: `δ_n-δ_{n+1}≥
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

★★★**結論(正直な評価)**: `/goal` の 13/13 は、このセッションでの
継続作業だけでは現実的な時間内に到達できない規模の作業(Theorem 1.2
1件でさらに大きな作業、§2-4 は新規ライブラリ構築に近い規模)である
可能性が高い。ただし Lemma 1.1 で確立した一般補題(especially
`subsingleton_H1Cotangent_self`・`polynomialKaehlerSplit`)は Theorem
1.2 の帰納法でも再利用できる見込みがあり、**無駄になっていない**。

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
